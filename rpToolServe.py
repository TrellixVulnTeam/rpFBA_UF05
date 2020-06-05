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

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

#logging.disable(logging.INFO)
#logging.disable(logging.WARNING)

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
#singleFBA_hdd
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
                  dont_merge=True,
                  pathway_id='rp_pathway',
                  objective_id=None,
                  compartment_id='MNXC3',
                  fill_orphan_species=False,
                  species_group_id='central_species'):
    rpsbml = rpSBML.rpSBML(file_name)
    rpsbml.readSBML(sbml_path)
    #Save the central species
    groups = rpsbml.model.getPlugin('groups')
    central = groups.getGroup(species_group_id)
    cent_spe = [str(i.getIdRef()) for i in central.getListOfMembers()]
    logging.info('old central species: '+str(cent_spe))
    #rpsbml_gem = rpSBML.rpSBML(file_name, libsbml.readSBMLFromString(gem_sbml))
    rpsbml_gem = rpSBML.rpSBML(file_name)
    rpsbml_gem.readSBML(gem_sbml)
    rpsbml.mergeModels(rpsbml_gem)
    #TO REMOVE
    #rpsbml_gem.modelName = 'test'
    #rpsbml_gem.writeSBML(tmpOutputFolder)
    #rpfba = rpFBA.rpFBA(rpsbml_gem)
    ####### fraction of reaction ######
    if sim_type=='fraction':
        rpfba.runFractionReaction(source_reaction, source_coefficient, target_reaction, target_coefficient, fraction_of, is_max, pathway_id, objective_id)
    ####### FBA ########
    elif sim_type=='fba':
        rpfba.runFBA(target_reaction, target_coefficient, is_max, pathway_id, objective_id)
    ####### pFBA #######
    elif sim_type=='pfba':
        rpfba.runParsimoniousFBA(target_reaction, target_coefficient, fraction_of, is_max, pathway_id, objective_id)
    else:
        logging.error('Cannot recognise sim_type: '+str(sim_type))
    '''
    ###### multi objective #####
    elif sim_type=='multi_fba':
        rpfba.runMultiObjective(reactions, coefficients, is_max, pathway_id)
    '''
    if dont_merge:
        logging.info('Returning model with heterologous pathway only')
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
            #logging.info('Replacing group id: '+str(source_group.getId()))
            if source_group.getId()==species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                #TODO: #### replace the new potentially incorect central species with the normal ones #####
                #delete all the previous members
                logging.info('Removing central_species')
                for i in range(target_group.getNumMembers()):
                    logging.info('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in cent_spe:
                    logging.info('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)  
            elif source_group.getId() in target_groupsID:
                target_group = target_groups.getGroup(source_group.getId())
                target_group.setAnnotation(source_group.getAnnotation())
        #### add objectives ####
        source_fbc = rpfba.rpsbml.model.getPlugin('fbc')
        target_fbc = rpsbml.model.getPlugin('fbc')
        target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
        for source_obj in source_fbc.getListOfObjectives():
            source_obj_id = source_obj.getId()
            if source_obj.getId() in target_objID:
                target_obj = target_fbc.getObjective(source_obj.getId())
                target_obj.setAnnotation(source_obj.getAnnotation())
                for target_fluxObj in target_obj.getListOfFluxObjectives():
                    for source_fluxObj in source_obj.getListOfFluxObjectives():
                        if target_fluxObj.getReaction()==source_fluxObj.getReaction():
                            target_fluxObj.setAnnotation(source_fluxObj.getAnnotation())
            else:
                target_fbc.addObjective(source_obj)
        #rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
        target_fbc.setActiveObjectiveId(source_obj_id) #tmp random assigenement of objective
        rpsbml.writeSBML(tmpOutputFolder)
    else:
        logging.info('Returning the full model')
        rpfba.rpsbml.writeSBML(tmpOutputFolder)


##
#
#
def runFBA_hdd(inputTar,
               gem_sbml,
               outputTar,
               sim_type,
               source_reaction,
               target_reaction,
               source_coefficient,
               target_coefficient,
               isMax,
               fraction_of,
               dont_merge=True,
               pathway_id='rp_pathway',
               objective_id=None,
               compartment_id='MNXC3',
               fill_orphan_species=False,
               species_group_id='central_species'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            #open the model as a string
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                logging.info('############## '+str(fileName)+' ################')
                try:
                    #logging.info('Running single FBA with the following parameters:')
                    #logging.info('\t')
                    singleFBA_hdd(fileName,
                                  sbml_path,
                                  gem_sbml,
                                  sim_type,
                                  source_reaction,
                                  target_reaction,
                                  source_coefficient,
                                  target_coefficient,
                                  isMax,
                                  fraction_of,
                                  tmpOutputFolder,
                                  dont_merge,
                                  pathway_id,
                                  objective_id,
                                  compartment_id,
                                  fill_orphan_species,
                                  species_group_id)
                except OSError as e:
                    logging.warning(e)
                    logging.warning('Segmentation fault by Cobrapy')
                    pass
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpFBA has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))+'.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True




##
#
#
def runFBA_multi(inputTar,
                 gem_sbml,
                 outputTar,
                 sim_type,
                 source_reaction,
                 target_reaction,
                 source_coefficient,
                 target_coefficient,
                 is_max,
                 fraction_of,
                 dont_merge=True,
                 num_workers=10,
                 pathway_id='rp_pathway',
                 objective_id=None,
                 compartment_id='MNXC3',
                 fill_orphan_species=False,
                 species_group_id='central_species'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
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
                                                                     dont_merge,
                                                                     pathway_id,
                                                                     objective_id,
                                                                     compartment_id,
                                                                     fill_orphan_species,
                                                                     species_group_id,)))
            output = [p.get() for p in results]
            pool.close()
            pool.join()
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpFBA has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))+'.sbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True




##
#
#
def main(input_path,
         gem_sbml,
         output_path,
         sim_type,
         source_reaction,
         target_reaction,
         source_coefficient,
         target_coefficient,
         is_max,
         fraction_of,
         dont_merge=True,
         num_workers=10,
         pathway_id='rp_pathway',
         objective_id=None,
         compartment_id='MNXC3',
         fill_orphan_species=None,
         species_group_id='central_species'):
    #outputTar_obj = io.BytesIO()
    if num_workers==1:
        runFBA_hdd(input_path,
                   gem_sbml,
                   output_path,
                   str(sim_type),
                   str(source_reaction),
                   str(target_reaction),
                   float(source_coefficient),
                   float(target_coefficient),
                   is_max,
                   float(fraction_of),
                   bool(dont_merge),
                   str(pathway_id),
                   objective_id,
                   str(compartment_id),
                   str(species_group_id))
        return True
    elif num_workers>1:
        runFBA_multi(input_path,
                     gem_sbml,
                     output_path,
                     str(sim_type),
                     str(source_reaction),
                     str(target_reaction),
                     float(source_coefficient),
                     float(target_coefficient),
                     is_max,
                     float(fraction_of),
                     bool(dont_merge),
                     int(num_workers),
                     str(pathway_id),
                     objective_id,
                     str(compartment_id),
                     fill_orphan_species,
                     str(species_group_id))
        return True
    else:
        logging.error('Cannot have 0 or less workers: '+str(num_workers))
        return False
