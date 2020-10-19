import cobra
from cobra.flux_analysis import pfba
import libsbml
import tempfile
import glob

import logging


#TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes
class rpFBA:
    """Class to simulate an rpsbml object using different FBA types and objective functions
    """
    def __init__(self, rpsbml):
        """Default constructor

        :param rpsbml: The rpSBML object

        :type rpsbml: rpSBML
        """
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Started instance of rpFBA')
        self.rpsbml = rpsbml
        #TODO enable FBC if not done so
        self.cobraModel = None
        #self._convertToCobra()


    ##########################################################
    ################# Private Functions ######################
    ##########################################################


    def _checklibSBML(self, value, message):
        """Private function that checks the libSBML calls.

        Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html

        :param value: The libSBML command returned int
        :param message: The string that describes the call

        :type value: int
        :type message: str

        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :return: None
        :rtype: None
        """
        if value is None:
            self.logger.error('LibSBML returned a null value trying to ' + message + '.')
            raise AttributeError
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                self.logger.error(err_msg)
                raise AttributeError
        else:
            return None


    def _convertToCobra(self):
        """Convert the rpSBML object to cobra object

        :return: Success or failure of the function
        :rtype: bool
        """
        try:
            with tempfile.TemporaryDirectory() as tmpOutputFolder:
                self.rpsbml.writeSBML(tmpOutputFolder)
                #logging.info(glob.glob(tmpOutputFolder+'/*'))
                #logging.info(cobra.io.validate_sbml_model(glob.glob(tmpOutputFolder+'/*')[0]))
                self.cobraModel = cobra.io.read_sbml_model(glob.glob(tmpOutputFolder+'/*')[0], use_fbc_package=True)
            #self.cobraModel = cobra.io.read_sbml_model(self.rpsbml.document.toXMLNode().toXMLString(), use_fbc_package=True)
            #use CPLEX
            # self.cobraModel.solver = 'cplex'
        except cobra.io.sbml.CobraSBMLError as e:
            self.logger.error(e)
            self.logger.error('Cannot convert the libSBML model to Cobra')
            return False
        return True


    ##########################################################
    ################# Helper functions #######################
    ##########################################################


    #TODO: move this to rpSBML
    def writeAnalysisResults(self, objective_id, cobra_results, pathway_id='rp_pathway'):
        """Method to harcode into BRSynth annotations the results of a COBRA analysis

        :param objective_id: The id of the objective to optimise
        :param cobra_results: The cobrapy results object 
        :param pathway_id: The id of the heterologous pathway group (Default: rp_pathway)
    
        :type cobra_results: cobra.ModelSummary
        :type objective_id: str
        :type pathway_id: str
    
        :return: None
        :rtype: None
        """
        self.logger.debug('----- Setting the results for '+str(objective_id)+ ' -----')
        groups = self.rpsbml.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting groups plugin')
        rp_pathway = groups.getGroup(pathway_id)
        if rp_pathway==None:
            self.logger.warning('The group '+str(pathway_id)+' does not exist... creating it')
            self.rpsbml.createPathway(pathway_id)
            rp_pathway = groups.getGroup(pathway_id)
        self._checklibSBML(rp_pathway, 'Getting RP pathway')
        #write the results to the rp_pathway
        self.logger.debug('Set '+str(pathway_id)+' with '+str('fba_'+str(objective_id))+' to '+str(cobra_results.objective_value))
        self.rpsbml.addUpdateBRSynth(rp_pathway, 'fba_'+str(objective_id), str(cobra_results.objective_value), 'mmol_per_gDW_per_hr', False)
        #get the objective
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC plugin')
        obj = fbc_plugin.getObjective(objective_id)
        self._checklibSBML(obj, 'Getting objective '+str(objective_id))
        self.rpsbml.addUpdateBRSynth(obj, 'flux_value', str(cobra_results.objective_value), 'mmol_per_gDW_per_hr', False)
        self.logger.debug('Set the objective '+str(objective_id)+' a flux_value of '+str(cobra_results.objective_value))
        for flux_obj in obj.getListOfFluxObjectives():
            #sometimes flux cannot be returned
            if cobra_results.fluxes.get(flux_obj.getReaction())==None:
                self.logger.warning('Cobra BUG: Cannot retreive '+str(flux_obj.getReaction())+' flux from cobrapy... setting to 0.0')
                self.rpsbml.addUpdateBRSynth(flux_obj, 'flux_value', str(0.0), 'mmol_per_gDW_per_hr', False)
                self.logger.debug('Set the reaction '+str(flux_obj.getReaction())+' a flux_value of '+str(0.0))
            else:
                self.rpsbml.addUpdateBRSynth(flux_obj, 'flux_value', str(cobra_results.fluxes.get(flux_obj.getReaction())), 'mmol_per_gDW_per_hr', False)
                self.logger.debug('Set the reaction '+str(flux_obj.getReaction())+' a flux_value of '+str(cobra_results.fluxes.get(flux_obj.getReaction())))
        #write all the results to the reactions of pathway_id
        for member in rp_pathway.getListOfMembers():
            reac = self.rpsbml.model.getReaction(member.getIdRef())
            if reac==None:
                self.logger.error('Cannot retreive the following reaction: '+str(member.getIdRef()))
                #return False
                continue
            self.logger.debug('Set the reaction '+str(member.getIdRef())+' a '+str('fba_'+str(objective_id))+' of '+str(cobra_results.fluxes.get(reac.getId())))
            self.rpsbml.addUpdateBRSynth(reac, 'fba_'+str(objective_id), str(cobra_results.fluxes.get(reac.getId())), 'mmol_per_gDW_per_hr', False)


    ##################################################################
    ######################## Model runs ##############################
    ##################################################################


    def runMultiObjective(self,
                          reactions,
                          coefficients,
                          is_max=True,
                          pathway_id='rp_pathway',
                          objective_id=None):
        """Run FBA using multiple objectives

        :param reactions: The ids of the reactions involved in the objective
        :param coefficients: The coefficients associated with the reactions id
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)
    
        :type reactions: list
        :type coefficients: list
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str
    
        :return: Success or failure of the function
        :rtype: bool
        """
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        objective_id = self.rpsbml.findCreateObjective(reactions, coefficients, is_max)
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self._convertToCobra():
            return False
        cobra_results = self.cobraModel.optimize()
        self.writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return True


    def runFBA(self, reaction_id, coefficient=1.0, is_max=True, pathway_id='rp_pathway', objective_id=None):
        """Run FBA using a single objective

        :param reaction_id: The id of the reactions involved in the objective
        :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)
    
        :type reaction_id: str
        :type coefficient: float
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str
    
        :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
        :rtype: tuple
        """
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        objective_id = self.rpsbml.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self._convertToCobra():
            return 0.0, False
        cobra_results = self.cobraModel.optimize()
        self.writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return cobra_results.objective_value, True


    def runParsimoniousFBA(self, reaction_id, coefficient=1.0, fraction_of_optimum=0.95, is_max=True, pathway_id='rp_pathway', objective_id=None):
        """Run parsimonious FBA using a single objective

        :param reaction_id: The id of the reactions involved in the objective
        :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
        :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.95)
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)
    
        :type reaction_id: str
        :type coefficient: float
        :type fraction_of_optimum: float
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str

        :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
        :rtype: tuple
        """
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        objective_id = self.rpsbml.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self._convertToCobra():
            return 0.0, False
        cobra_results = pfba(self.cobraModel, fraction_of_optimum)
        self.writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return cobra_results.objective_value, True


    def runFractionReaction(self, 
                            source_reaction, 
                            source_coefficient, 
                            target_reaction, 
                            target_coefficient, 
                            fraction_of_source=0.75, 
                            is_max=True, 
                            pathway_id='rp_pathway', 
                            objective_id=None):
        """Optimise for a target reaction while fixing a source reaction to the fraction of its optimum

        :param source_reaction: The id of the source reaction
        :param source_coefficient: The source coefficient associated with the source reaction id
        :param target_reaction: The id of the target reaction
        :param target_coefficient: The source coefficient associated with the target reaction id
        :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.75)
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)
    
        :type source_reaction: str
        :type source_coefficient: float
        :type target_reaction: str
        :type target_coefficient: float
        :type fraction_of_optimum: float
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str

        :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
        :rtype: tuple
        """
        #retreive the biomass objective and flux results and set as maxima
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        self.logger.debug('findCreateObjective: '+str(source_reaction))
        source_obj_id = self.rpsbml.findCreateObjective([source_reaction], [source_coefficient], is_max)
        #TODO: use the rpSBML BRSynth annotation parser
        source_flux = None
        try:
            fbc_obj = fbc_plugin.getObjective(source_obj_id)
            #TODO: if this is None need to set it up 
            fbc_obj_annot = fbc_obj.getAnnotation()
            if not fbc_obj_annot:
                raise ValueError
            source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
            self.logger.debug('Already calculated flux for '+str(source_obj_id))
        except (AttributeError, ValueError) as e:
            self.logger.debug('Performing FBA to calculate the source reaction')
            ### FBA ###
            #self.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
            self._checklibSBML(fbc_plugin.setActiveObjectiveId(source_obj_id),
                    'Setting active objective '+str(source_obj_id))
            if not self._convertToCobra():
                self.logger.error('Converting libSBML to CobraPy returned False')
                self.writeAnalysisResults(source_obj_id, 0.0, pathway_id)
                return 0.0, False
            cobra_results = self.cobraModel.optimize()
            self.writeAnalysisResults(source_obj_id, cobra_results, pathway_id)
            # cobra_results.objective_value
            fbc_obj = fbc_plugin.getObjective(source_obj_id)
            fbc_obj_annot = fbc_obj.getAnnotation()
            if fbc_obj_annot==None:
                self.logger.error('No annotation available for: '+str(source_obj_id))
            source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
        #TODO: add another to check if the objective id exists
        self.logger.debug('FBA source flux ('+str(source_reaction)+') is: '+str(source_flux))
        if not objective_id:
            objective_id = 'obj_'+str(target_reaction)+'__restricted_'+str(source_reaction)
        #self.logger.debug('findCreateObjective() for '+str(objective_id))
        objective_id = self.rpsbml.findCreateObjective([target_reaction], [target_coefficient], is_max, objective_id)
        self.logger.debug('Optimising the objective: '+str(objective_id))
        self.logger.debug('Setting upper bound: '+str(source_flux*fraction_of_source))
        self.logger.debug('Setting loer bound: '+str(source_flux*fraction_of_source))
        old_upper_bound, old_lower_bound = self.rpsbml.setReactionConstraints(source_reaction,
                                                                              source_flux*fraction_of_source,
                                                                              source_flux*fraction_of_source)
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self._convertToCobra():
            self.logger.error('Converting libSBML to CobraPy returned False')
            #although this may not be the greatest idea, set flux to 0.0 when cobrapy error
            self.writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobraModel.optimize()
        self.writeAnalysisResults(objective_id, cobra_results, pathway_id)
        ##### print the biomass results ######
        #self.logger.debug('Biomass: '+str(cobra_results.fluxes.biomass))
        #self.logger.debug('Target: '+str(cobra_results.fluxes.RP1_sink))
        #reset the bounds to the original values for the target
        old_upper_bound, old_lower_bound = self.rpsbml.setReactionConstraints(source_reaction,
                                                                              old_upper_bound,
                                                                              old_lower_bound)
        self.logger.debug('The objective '+str(objective_id)+' results '+str(cobra_results.objective_value))
        return cobra_results.objective_value, True




    ########################################################################
    ############################### FBA pathway ranking ####################
    ########################################################################

    #1) Number of interventions
    # need to calculate the number of steps that are not native to know the number of interventions

    #2) Maximal growth rate

    #3) Minimum product yeild at maximal growth rate

    #4) Minimum product yeild

    #5) Anaerobic condition

    #6) Number of potentially disruptive products

        #Toxicity?

    #7) Number of accessible metabolites (avoid intermediate accumulation)

    #8) Thermodynamics (MDF)

    #9) The overlap of the same changes --> might not be applicable in our case

    #10) Reduced model

    #11) ECM
