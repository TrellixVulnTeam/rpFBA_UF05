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

import json
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api
import tempfile

sys.path.insert(0, '/home/')
import rpTool as rpFBA

import rpSBML

###################################################################################
###################################################################################
###################################################################################

import inspect
import traceback

from functools import wraps
from multiprocessing import Process, Queue


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

#hack to stop the memory leak. Indeed it seems that looping through rpFBA and the rest causes a memory leak... According to: https://github.com/opencobra/cobrapy/issues/568 there is still memory leak issues with cobrapy. looping through hundreds of models and running FBA may be the culprit

########################## use RAM ######################
"""
'''
##
#
#
@processify
def singleFBA_mem(fileName,
                  sbml_path,
                  inModel_string,
                  source_reaction,
                  target_reaction,
                  fraction_of_source,
                  tmpOutputFolder,
                  dontMerge,
                  pathway_id,
                  fill_orphan_species,
                  compartment_id):
    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(rpsbml_string))
    input_rpsbml = rpSBML.rpSBML(fileName, libsbml.readSBMLFromString(inModel_string))
    rpsbml.mergeModels(input_rpsbml, pathway_id, fill_orphan_species, compartment_id)
    rpfba = rpFBA.rpFBA(input_rpsbml)
    #rpfba.allObj(pathway_id)
    rpfba.runFractionReaction(source_reaction, target_reaction, fraction_of_source, pathway_id)
    if dontMerge:
        ##### pass FBA results to the original model ####
        groups = rpfba.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        for member in rp_pathway.getListOfMembers():
            #### reactions
            reacFBA = rpfba.rpsbml.model.getReaction(member.getIdRef())
            reacIN = rpsbml.model.getReaction(member.getIdRef())
            reacIN.setAnnotation(reacFBA.getAnnotation())
            #### species TODO: only for shadow price
        #### species #TODO: implement this for shadow prices
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
        return libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
    else:
        return libsbml.writeSBMLToString(input_rpsbml.document).encode('utf-8')


##
#
#
def runFBA_mem(inputTar, inModel_bytes, outTar, dontMerge, pathway_id='rp_pathway', fill_orphan_species=False, compartment_id='MNXC3'):
    #open the model as a string
    inModel_string = inModel_bytes.read().decode('utf-8')
    #loop through all of them and run FBA on them
    with tarfile.open(fileobj=outTar, mode='w:xz') as tf:
        with tarfile.open(fileobj=inputTar, mode='r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    data = singleFBA_mem(member.name, in_tf.extractfile(member).read().decode("utf-8"), inModel_string, dontMerge, pathway_id)
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)
'''
"""

####################### use HDD ############################

##
#
#
@processify
def singleFBA_hdd(fileName,
                  sbml_path,
                  inModel_string,
                  sim_type,
                  reactions,
                  coefficients,
                  isMax,
                  fraction_of,
                  tmpOutputFolder,
                  dontMerge,
                  pathway_id,
                  fill_orphan_species,
                  compartment_id):
    rpsbml = rpSBML.rpSBML(fileName)
    rpsbml.readSBML(sbml_path)
    input_rpsbml = rpSBML.rpSBML(fileName, libsbml.readSBMLFromString(inModel_string))
    rpsbml.mergeModels(input_rpsbml, pathway_id, fill_orphan_species, compartment_id)
    rpfba = rpFBA.rpFBA(input_rpsbml)
    ####### fraction of reaction ######
    if sim_type=='fraction':
        rpfba.runFractionReaction(reactions[0], reactions[1], fraction_of, pathway_id)
    ####### FBA ########
    elif sim_type=='fba':
        rpfba.runFBA(reactions[0], isMax, pathway_id)
    ####### pFBA #######
    elif sim_type=='pfba':
        rpfba.runParsimoniousFBA(reactions[0], fraction_of, isMax, pathway_id)
    ###### multi objective #####
    elif sim_type=='multi_fba':
        rpfba.runMultiObjective(reactions, coefficients, isMax, pathway_id)
    else:
        logging.error('Cannot recognise sim_type: '+str(sim_type))    
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
               reactions,
               coefficients,
               isMax,
               fraction_of,
               dontMerge,
               pathway_id='rp_pathway',
               fill_orphan_species=False,
               compartment_id='MNXC3'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            #open the model as a string
            inModel_string = inModel_bytes.read().decode('utf-8')
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                singleFBA_hdd(fileName,
                              sbml_path,
                              inModel_string,
                              sim_type,
                              reactions,
                              coefficients,
                              isMax,
                              float(fraction_of),
                              tmpOutputFolder,
                              dontMerge,
                              pathway_id,
                              fill_orphan_species,
                              compartment_id)
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', ''))+'.rpsbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))


#######################################################
############## REST ###################################
#######################################################


app = Flask(__name__)
api = Api(app)


def stamp(data, status=1):
    appinfo = {'app': 'rpFBA', 'version': '1.0',
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(),
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp(Resource):
    """ REST App."""
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):
        inputTar = request.files['inputTar']
        inSBML = request.files['inSBML']
        params = json.load(request.files['data'])
        #pass the files to the rpReader
        outputTar = io.BytesIO()
        #### MEM ####
        '''
        runFBA_mem(inputTar,
                   inSBML,
                   outputTar,
                   str(params['source_reaction']),
                   str(params['target_reaction']),
                   float(params['fraction_of_source'])
                   bool(params['dontMerge']),
                   str(params['pathway_id']),
                   bool(params['fill_orphan_species']),
                   str(params['compartment_id']))
        '''
        #### HDD ####
        runFBA_hdd(inputTar,
                   inSBML,
                   outputTar,
                   str(params['sim_type']),
                   list(params['reactions']),
                   list(params['coefficients']),
                   bool(params['isMax']),
                   float(params['fraction_of']),
                   bool(params['dontMerge']),
                   str(params['pathway_id']),
                   bool(params['fill_orphan_species']),
                   str(params['compartment_id']))
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpFBA.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
