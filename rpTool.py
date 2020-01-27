import cobra
from cobra.flux_analysis import pfba
import libsbml

import logging

## Class to simulate an rpsbml object using different FBA types and objective functions
#
# At this point we want to have the BIOMASS, target and shared objective
#TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes
class rpFBA:
    def __init__(self, rpsbml):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpFBA')
        self.rpsbml = rpsbml
        #TODO enable FBC if not done so
        self.cobraModel = None
        self._convertToCobra()


    ##########################################################
    ################# Private Functions ######################
    ##########################################################


    ## Check the libSBML calls
    #
    # Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
    #
    # @param value The SBML call
    # @param message The string that describes the call
    def _checklibSBML(self, value, message):
        if value is None:
            self.logger.error('LibSBML returned a null value trying to ' + message + '.')
            raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                self.logger.error(err_msg)
                raise SystemExit(err_msg)
        else:
            #self.logger.info(message)
            return None


    ## Pass the libSBML file to Cobra
    #
    #
    def _convertToCobra(self):
        try:
            self.cobraModel = cobra.io.read_sbml_model(self.rpsbml.document.toXMLNode().toXMLString(),
                    use_fbc_package=True)
            #use CPLEX
            # self.cobraModel.solver = 'cplex'
        except cobra.io.sbml.CobraSBMLError as e:
            self.logger.error(e)
            self.logger.error('Cannot convert the libSBML model to Cobra')


    ##########################################################
    ################# Helper functions #######################
    ##########################################################


    ## Method to harcode into BRSynth annotations the results of a COBRA analysis
    #
    #
    def writeAnalysisResults(self, fbc_obj, cobra_results, pathway_id='rp_pathway'):
        groups = self.rpsbml.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting groups plugin')
        rp_pathway = groups.getGroup(pathway_id)
        self._checklibSBML(rp_pathway, 'Getting RP pathway')
        self.rpsbml.addUpdateBRSynth(fbc_obj, 'flux_value', str(cobra_results.objective_value), 'mmol_per_gDW_per_hr', False)
        for flux_obj in fbc_obj.getListOfFluxObjectives():
            self.rpsbml.addUpdateBRSynth(flux_obj, 'flux_value', str(cobra_results.objective_value), 'mmol_per_gDW_per_hr', False)
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                if reac==None:
                    self.logger.error('Cannot retreive the following reaction: '+str(member.getIdRef()))
                    return False
                self.rpsbml.addUpdateBRSynth(reac, 'fba_'+str(fbc_obj.getId()), str(cobra_results.fluxes.get(reac.getId())), 'mmol_per_gDW_per_hr', False)


    ## set Bi-objective 
    #
    #
    def setMultiObjective(self,
                          reactions,
                          coefficients=[0.5, 0.5],
                          fluxobj_id='rpFBA_biObj',
                          isMax=True):
        try:
            self.cobraModel.get_by_any(biomass_reaction)
        except KeyError:
            self.logger.error('Cannot find the reaction '+str(biomass_reaction))
            return False
        self.model.createMultiFluxObj(fluxobj_id,
                                      [rpFBA_sink_reaction, biomass_reaction],
                                      coefficients,
                                      isMax)


    ##################################################################
    ######################## Model runs ##############################
    ##################################################################


    ##
    #
    #
    def runAllParsimoniousFBA(self, fraction_of_optimum=0.95, pathway_id='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        for obj in fbc_plugin.getListOfObjectives():
            self.runParsimoniousFBA(obj.getId(), fraction_of_optimum, pathway_id)


    ## 
    #
    #
    def runAllFBA(self, pathway_id='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        for obj in fbc_plugin.getListOfObjectives():
            self.runFBA(obj.getId(), pathway_id)


    ##
    #
    #
    def runFBA(self, objective_id, pathway_id='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        fbc_obj = fbc_plugin.getObjective(objective_id)
        self._checklibSBML(fbc_obj, 'Getting the following objective '+str(objective_id))
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        self._convertToCobra()
        cobra_results = self.cobraModel.optimize()
        self.writeAnalysisResults(fbc_obj, cobra_results, pathway_id)


    ##
    #
    #
    def runParsimoniousFBA(self, objective_id, fraction_of_optimum=0.95, pathway_id='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        fbc_obj = fbc_plugin.getObjective(objective_id)
        self._checklibSBML(fbc_obj, 'Getting the following objective '+str(objective_id))
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        self._convertToCobra()
        cobra_results = pfba(self.cobraModel, fraction_of_optimum)
        writeAnalysisResults(fbc_obj, cobra_results, pathway_id)


    ## Optimise for a target reaction while fixing a source reaction to the fraction of its optimum
    #
    #
    def runFractionReaction(self, source_reaction, target_reaction, fraction_of_source=0.75, pathway_id='rp_pathway'):
        #retreive the biomass objective and flux results and set as maxima
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        source_obj_id = self.rpsbml.findCreateObjective(source_reaction, pathway_id)
        #target_obj_id = self.rpsbml.findCreateObjective(target_reaction, pathway_id)
        source_flux = None
        fbc_obj = fbc_plugin.getObjective(source_obj_id)
        fbc_obj_annot = fbc_obj.getAnnotation()
        if not fbc_obj_annot==None:
            if not fbc_obj_annot.getAnnotationString()=='':
                source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
            else:
                self.runFBA(source_obj_id, pathway_id)
                fbc_obj = fbc_plugin.getObjective(source_obj_id)
                fbc_obj_annot = fbc_obj.getAnnotation()
                if fbc_obj_annot==None:
                    self.logger.error('There is an error getting the flux for source objective: '+str(source_reaction))
                    return 0.0
                source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
        else:
            self.runFBA(source_obj_id, pathway_id)
            fbc_obj = fbc_plugin.getObjective(source_obj_id)
            fbc_obj_annot = fbc_obj.getAnnotation()
            if fbc_obj_annot==None:
                self.logger.error('There is an error getting the flux for source objective: '+str(source_reaction))
                return 0.0
            source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
        #set bounds biomass as a fraction
        target_obj_id = str(target_reaction)+'__restricted_'+str(source_reaction)
        self.rpsbml.createMultiFluxObj(str(target_reaction)+'__restricted_'+str(source_reaction), ['RP1_sink'], [1])
        old_upper_bound, old_lower_bound = self.rpsbml.setReactionConstraints(source_reaction,
                                                                              source_flux*fraction_of_source,
                                                                              source_flux*fraction_of_source)
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(target_obj_id),
                'Setting active objective '+str(target_obj_id))
        self._convertToCobra()
        fbc_obj = fbc_plugin.getObjective(target_obj_id)
        self._checklibSBML(fbc_obj, 'Getting the following objective '+str(target_obj_id))
        cobra_results = self.cobraModel.optimize()
        self.writeAnalysisResults(fbc_obj, cobra_results)
        #reset the bounds to the original values for the target
        old_upper_bound, old_lower_bound = self.rpsbml.setReactionConstraints(source_reaction,
                                                                              old_upper_bound,
                                                                              old_lower_bound)
        return cobra_results.objective_value




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
