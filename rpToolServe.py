#!/usr/bin/env python3

#from contextlib import closing
#import time
import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import io
import tarfile
import glob
import logging
import tempfile
import shutil

sys.path.insert(0, '/home/')
import rpTool as rpFBA

import rpSBML


###################################################################################
################################## processify #####################################
###################################################################################

#hack to stop the memory leak. Indeed it seems that looping through rpFBA and the rest causes a memory leak... According to: https://github.com/opencobra/cobrapy/issues/568 there is still memory leak issues with cobrapy. looping through hundreds of models and running FBA may be the culprit

import inspect
import traceback
import signal
from functools import wraps
from multiprocessing import Process, Queue

'''
This is to deal with an error caused by Cobrapy segmentation fault
'''
def handler(signum, frame):
    raise OSError('CobraPy is throwing a segmentation fault')


class Sentinel:
    pass


def processify(func):
    '''Decorator to run a function as a process.
    Be sure that every argument and the return value
    is *pickable*.
    The created process is joined, so the code does not
    run in parallel.
    '''

    def process_generator_func(q, *args, **kwargs):
        result = None
        error = None
        it = iter(func())
        while error is None and result != Sentinel:
            try:
                result = next(it)
                error = None
            except StopIteration:
                result = Sentinel
                error = None
            except Exception:
                ex_type, ex_value, tb = sys.exc_info()
                error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
                result = None
            q.put((result, error))

    def process_func(q, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
        except Exception:
            ex_type, ex_value, tb = sys.exc_info()
            error = ex_type, ex_value, ''.join(traceback.format_tb(tb))
            result = None
        else:
            error = None

        q.put((result, error))

    def wrap_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_func.__name__ = func.__name__ + 'processify_func'
        setattr(sys.modules[__name__], process_func.__name__, process_func)

        signal.signal(signal.SIGCHLD, handler) #This is to catch the segmentation error

        q = Queue()
        p = Process(target=process_func, args=[q] + list(args), kwargs=kwargs)
        p.start()
        result, error = q.get()
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

        return result

    def wrap_generator_func(*args, **kwargs):
        # register original function with different name
        # in sys.modules so it is pickable
        process_generator_func.__name__ = func.__name__ + 'processify_generator_func'
        setattr(sys.modules[__name__], process_generator_func.__name__, process_generator_func)

        signal.signal(signal.SIGCHLD, handler) #This is to catch the segmentation error

        q = Queue()
        p = Process(target=process_generator_func, args=[q] + list(args), kwargs=kwargs)
        p.start()

        result = None
        error = None
        while error is None:
            result, error = q.get()
            if result == Sentinel:
                break
            yield result
        p.join()

        if error:
            ex_type, ex_value, tb_str = error
            message = '%s (in subprocess)\n%s' % (str(ex_value), tb_str)
            raise ex_type(message)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if inspect.isgeneratorfunction(func):
            return wrap_generator_func(*args, **kwargs)
        else:
            return wrap_func(*args, **kwargs)
    return wrapper


###########################################################
################## multiprocesses run #####################
###########################################################

#This is a non-deamonic multiprocessing method that can be used in combination with processify

import multiprocessing
import multiprocessing.pool
import time


class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class nonDeamonicPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
    


####################### use HDD ############################


##
#
#
@processify
def singleFBA_hdd(file_name,
                  sbml_path,
                  gem_sbml,
                  sim_type,
                  source_reaction,
                  target_reaction,
                  source_coefficient,
                  target_coefficient,
                  is_max,
                  fraction_of,
                  tmpOutputFolder,
                  dontMerge,
                  pathway_id,
                  compartment_id,
                  fill_orphan_species):
    logging.warning(sbml_path)
    logging.warning(file_name)
    logging.warning('######################')
    rpsbml = rpSBML.rpSBML(file_name)
    rpsbml.readSBML(sbml_path)
    #input_rpsbml = rpSBML.rpSBML(file_name, libsbml.readSBMLFromString(gem_sbml))
    input_rpsbml = rpSBML.rpSBML(file_name)
    input_rpsbml.readSBML(gem_sbml)
    rpsbml.mergeModels(input_rpsbml)
    rpfba = rpFBA.rpFBA(input_rpsbml)
    ####### fraction of reaction ######
    if sim_type=='fraction':
        rpfba.runFractionReaction(source_reaction, target_reaction, fraction_of, pathway_id)
    ####### FBA ########
    elif sim_type=='fba':
        rpfba.runFBA(target_reaction, is_max, pathway_id)
    ####### pFBA #######
    elif sim_type=='pfba':
        rpfba.runParsimoniousFBA(target_reaction, fraction_of, is_max, pathway_id)
    else:
        logging.error('Cannot recognise sim_type: '+str(sim_type))    
    '''
    ###### multi objective #####
    elif sim_type=='multi_fba':
        rpfba.runMultiObjective(reactions, coefficients, is_max, pathway_id)
    '''
    if dontMerge:
        groups = rpfba.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        for member in rp_pathway.getListOfMembers():
            #### reaction annotation
            reacFBA = rpfba.rpsbml.model.getReaction(member.getIdRef())
            reacIN = rpsbml.model.getReaction(member.getIdRef())
            reacIN.setAnnotation(reacFBA.getAnnotation())
            #### species TODO: only for shadow price
        #### add groups ####
        source_groups = rpfba.rpsbml.model.getPlugin('groups')
        target_groups = rpsbml.model.getPlugin('groups')
        target_groupsID = [i.getId() for i in target_groups.getListOfGroups()]
        for source_group in source_groups.getListOfGroups():
            if source_group.getId() in target_groupsID:
                target_group = target_groups.getGroup(source_group.getId())
                target_group.setAnnotation(source_group.getAnnotation())
        #### add objectives ####
        source_fbc = rpfba.rpsbml.model.getPlugin('fbc')
        target_fbc = rpsbml.model.getPlugin('fbc')
        target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
        for source_obj in source_fbc.getListOfObjectives():
            if source_obj.getId() in target_objID:
                target_obj = target_fbc.getObjective(source_obj.getId())
                target_obj.setAnnotation(source_obj.getAnnotation())
                for target_fluxObj in target_obj.getListOfFluxObjectives():
                    for source_fluxObj in source_obj.getListOfFluxObjectives():
                        if target_fluxObj.getReaction()==source_fluxObj.getReaction():
                            target_fluxObj.setAnnotation(source_fluxObj.getAnnotation())
            else:
                target_fbc.addObjective(source_obj)
        rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
        rpsbml.writeSBML(tmpOutputFolder)
    else:
        rpfba.rpsbml.writeSBML(tmpOutputFolder)


##
#
#
def runFBA_hdd(inputTar,
               inModel_bytes,
               outputTar,
               sim_type,
               source_reaction,
               target_reaction,
               source_coefficient,
               target_coefficient,
               isMax,
               fraction_of,
               dontMerge,
               pathway_id='rp_pathway',
               compartment_id='MNXC3',
               fill_orphan_species=False):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            #open the model as a string
            inModel_string = inModel_bytes.read().decode('utf-8')
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                try:
                    singleFBA_hdd(fileName,
                                  sbml_path,
                                  inModel_string,
                                  sim_type,
                                  source_reaction,
                                  target_reaction,
                                  source_coefficient,
                                  target_coefficient,
                                  isMax,
                                  fraction_of,
                                  tmpOutputFolder,
                                  dontMerge,
                                  pathway_id,
                                  compartment_id,
                                  fill_orphan_species)
                except OSError as e:
                    logging.warning(e)
                    logging.warning('Segmentation fault by Cobrapy') 
                    pass
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))+'.rpsbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))




##
#
#
def runFBA_multi(inputTar,
                 inModel_bytes,
                 outputTar,
                 sim_type,
                 source_reaction,
                 target_reaction,
                 source_coefficient,
                 target_coefficient,
                 is_max,
                 fraction_of,
                 dontMerge,
                 num_workers=10,
                 pathway_id='rp_pathway',
                 compartment_id='MNXC3',
                 fill_orphan_species=False):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            #open the model as a string
            gem_sbml = inModel_bytes.read().decode('utf-8')
            #HERE SPECIFY THE NUMBER OF CORES
            pool = nonDeamonicPool(processes=num_workers)
            results = []
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                file_name = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                results.append(pool.apply_async(singleFBA_hdd, args=(file_name,
                                                                     sbml_path,
                                                                     gem_sbml,
                                                                     sim_type,
                                                                     source_reaction,
                                                                     target_reaction,
                                                                     source_coefficient,
                                                                     target_coefficient,
                                                                     is_max,
                                                                     fraction_of,
                                                                     tmpOutputFolder,
                                                                     dontMerge,
                                                                     pathway_id,
                                                                     compartment_id,
                                                                     fill_orphan_species,)))
            output = [p.get() for p in results]
            logging.info(output)
            pool.close()
            pool.join()
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))+'.rpsbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))




##
#
#
def main(input_path, 
         full_sbml_path, 
         output_path, 
         sim_type, 
         source_reaction, 
         target_reaction, 
         source_coefficient, 
         target_coefficient, 
         is_max, 
         fraction_of, 
         dont_merge, 
         num_workers,
         pathway_id, 
         compartment_id):
    with open(input_path, 'rb') as input_bytes:
        with open(full_sbml_path, 'rb') as full_sbml_bytes:
            outputTar_obj = io.BytesIO()
            '''
            runFBA_multi(input_bytes,
                         full_sbml_bytes,
                         outputTar_obj,
                         str(sim_type),
                         str(source_reaction),
                         str(target_reaction),
                         float(source_coefficient),
                         float(target_coefficient),
                         bool(is_max),
                         float(fraction_of),
                         bool(dont_merge),
                         int(num_workers),
                         str(pathway_id),
                         str(compartment_id))
            '''
            '''DEPRECATED
            '''
            runFBA_hdd(input_bytes,
                       full_sbml_bytes,
                       outputTar_obj,
                       str(sim_type),
                       str(source_reaction),
                       str(target_reaction),
                       float(source_coefficient),
                       float(target_coefficient),
                       bool(is_max),
                       float(fraction_of),
                       bool(dont_merge),
                       str(pathway_id),
                       str(compartment_id))
            '''
            '''
            ########## IMPORTANT #####
            outputTar_obj.seek(0)
            ##########################
            with open(output_path, 'wb') as f:
                shutil.copyfileobj(outputTar_obj, f, length=131072)

