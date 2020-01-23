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


    def convertToCobra(self):
        try:
            self.cobraModel = cobra.io.read_sbml_model(self.rpsbml.document.toXMLNode().toXMLString(),
                    use_fbc_package=True)
            #use CPLEX
            # self.cobraModel.solver = 'cplex'
        except cobra.io.sbml.CobraSBMLError as e:
            self.logger.error(e)
            self.logger.error('Cannot convert the libSBML model to Cobra')


    ## TODO: use the objective from the original model (GEM) that contains the biomass function
    # and return the flux for the biomass reaction. This value will be used to normalise the FBA
    # score
    #   
    def allObj(self, pathway_id='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        groups = self.rpsbml.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting groups plugin')
        rp_pathway = groups.getGroup(pathway_id)
        self._checklibSBML(rp_pathway, 'Getting RP pathway')
        #for objId in [i.getId() for i in fbc_plugin.getListOfObjectives()]:
        for fbc_obj in fbc_plugin.getListOfObjectives():
            #find all the species of the reactions in the rp_pathway to return the shadow price
            ''' TO BE DETERMINED IF USED 
            mem = []
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                for pro in reac.getListOfProducts():
                    mem.append(pro.getSpecies())
                for rea in reac.getListOfReactants():
                    mem.append(rea.getSpecies())
            '''
            #run the FBA
            self._checklibSBML(fbc_plugin.setActiveObjectiveId(fbc_obj.getId()), 'Setting active objective '+str(fbc_obj.getId()))
            self.convertToCobra()
            res = self.cobraModel.optimize()
            #add the objective value to the fbc obj
            tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:objective_value units="mmol_per_gDW_per_hr" value="'+str(res.objective_value)+'" /> </brsynth:brsynth>')
            try:
                obj_annot = fbc_obj.getAnnotation()
                brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                brsynth_annot.addChild(tmpAnnot.getChild('objective_value'))
            except AttributeError:
                annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:BRSynth rdf:about="#">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''
                fbc_obj.setAnnotation(annotation)
                obj_annot = fbc_obj.getAnnotation()
                brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                brsynth_annot.addChild(tmpAnnot.getChild('objective_value'))
            #add the results of FBA run to the annotation of FBA objective
            for flux_obj in fbc_obj.getListOfFluxObjectives():
                try:
                    obj_annot = flux_obj.getAnnotation()
                    brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:reac_obj_value units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(flux_obj.getReaction()))+'" /> </brsynth:brsynth>')
                    brsynth_annot.addChild(tmpAnnot.getChild('reac_obj_value'))
                except AttributeError:
                    annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
      xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
        <rdf:BRSynth rdf:about="#">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
                    flux_obj.setAnnotation(annotation)
                    obj_annot = flux_obj.getAnnotation()
                    brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:reac_obj_value units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(flux_obj.getReaction()))+'" /> </brsynth:brsynth>')
                    brsynth_annot.addChild(tmpAnnot.getChild('reac_obj_value'))
            #add the results of the FBA reactions for each rp_pathway reactions 
            #TODO: need to take care of the case where the annotation exists already --> overwrite
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                reac_annot = reac.getAnnotation()
                brsynth_annot = reac_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:fba_'+str(fbc_obj.getId())+' units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(member.getIdRef()))+'" /> </brsynth:brsynth>')
                brsynth_annot.addChild(tmpAnnot.getChild('fba_'+str(fbc_obj.getId())))
            ''' TO BE DETERMINED IF USED
            #update the shadow prices for species
            for speName in list(set(mem)): #remoce duplicates
                #only choose the heterologous species
                if len([x for x in speName.split('_') if x])==4:
                    spe = self.rpsbml.model.getSpecies(speName)
                    spe_annot = spe.getAnnotation()
                    brsynth_annot = spe_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:shadow_price_'+str(objId)+' units="mmol_per_gDW_per_hr" value="'+str(res.shadow_prices.get(speName))+'" /> </brsynth:brsynth>')
                    brsynth_annot.addChild(tmpAnnot.getChild('shadow_price_'+str(objId)))
            '''



    ########################################################################


    def writeAnalysisResults(self, fbc_obj, cobra_results, pathway_id='rp_pathway'):
        groups = self.rpsbml.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting groups plugin')
        rp_pathway = groups.getGroup(pathway_id)
        self._checklibSBML(rp_pathway, 'Getting RP pathway')
        #add the objective value to the fbc obj
        tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:flux_value units="mmol_per_gDW_per_hr" value="'+str(cobra_results.objective_value)+'" /> </brsynth:brsynth>')
        try:
            obj_annot = fbc_obj.getAnnotation()
            brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
            if brsynth_annot.getChild('flux_value').toXMLString()=='':
                brsynth_annot.addChild(tmpAnnot.getChild('flux_value'))
            else:
                result_brsynth_annot = brsynth_annot.getChild('flux_value')
                result_brsynth_annot.removeAttr('value')
                result_brsynth_annot.addAttr('value', str(cobra_results.objective_value))
        except AttributeError:
            annotation = '''<annotation>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
<rdf:BRSynth rdf:about="#">
  <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
  </brsynth:brsynth>
</rdf:BRSynth>
</rdf:RDF>
</annotation>'''
            fbc_obj.setAnnotation(annotation)
            obj_annot = fbc_obj.getAnnotation()
            brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
            brsynth_annot.addChild(tmpAnnot.getChild('flux_value'))
        #add the results of FBA run to the annotation of FBA objective
        for flux_obj in fbc_obj.getListOfFluxObjectives():
            try:
                obj_annot = flux_obj.getAnnotation()
                brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:flux_value units="mmol_per_gDW_per_hr" value="'+str(cobra_results.fluxes.get(flux_obj.getReaction()))+'" /> </brsynth:brsynth>')
            if brsynth_annot.getChild('flux_value').toXMLString()=='':
                brsynth_annot.addChild(tmpAnnot.getChild('flux_value'))
            else:
                result_brsynth_annot = brsynth_annot.getChild('flux_value')
                result_brsynth_annot.removeAttr('value')
                result_brsynth_annot.addAttr('value', str(cobra_results.fluxes.get(flux_obj.getReaction())))
            except AttributeError:
                annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:BRSynth rdf:about="#">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''
                flux_obj.setAnnotation(annotation)
                obj_annot = flux_obj.getAnnotation()
                brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:flux_value units="mmol_per_gDW_per_hr" value="'+str(cobra_results.fluxes.get(flux_obj.getReaction()))+'" /> </brsynth:brsynth>')
                brsynth_annot.addChild(tmpAnnot.getChild('flux_value'))
        #add the results of the FBA reactions for each rp_pathway reactions 
        #TODO: need to take care of the case where the annotation exists already --> overwrite
        for member in rp_pathway.getListOfMembers():
            reac = self.rpsbml.model.getReaction(member.getIdRef())
            reac_annot = reac.getAnnotation()
            brsynth_annot = reac_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
            tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:fba_'+str(fbc_obj.getId())+' units="mmol_per_gDW_per_hr" value="'+str(cobra_results.fluxes.get(member.getIdRef()))+'" /> </brsynth:brsynth>')
            if brsynth_annot.getChild('fba_'+str(fbc_obj.getId())).toXMLString()=='':
                brsynth_annot.addChild(tmpAnnot.getChild('fba_'+str(fbc_obj.getId())))
            else:
                result_brsynth_annot = brsynth_annot.getChild('fba_'+str(fbc_obj.getId()))
                result_brsynth_annot.removeAttr('value')
                result_brsynth_annot.addAttr('value', str(cobra_results.fluxes.get(flux_obj.getReaction())))




    
    ## Find the objective based on the reaction ID and if not found create and simulate it 
    #
    #
    def findCreateObjective(self, reaction_id, pathway_id='rp_pathway'):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        for objective in fbc_plugin.getListOfObjectives():
            if len(objective.getListOfFluxObjectives())==1:
                #try to retreive the results from the annotation
                flux_objective = objective.getListOfFluxObjectives()[0]
                if flux_objective.getReaction()==reaction_id:
                    try:
                        flux = float(flux_objective.getAnnotation().getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
                        return objective.getId(), float(flux)
                    except ValueError:
                        #need to run the flux for the biomass objective
                        self.runFBA(objective.getId(), pathway_id)
                        try:
                            flux = float(flux_objective.getAnnotation().getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
                            return objective.getId(), float(flux)
                        except ValueError:
                            self.logger.error('There is an error setting the annotation for source objective')
                            return '', 0.0
        #if here means that was not found -- create and run it
        self.rpsbml.createMultiFluxObj(str(reaction_id)+'_obj', [reaction_id], [1], True)
        self.runFBA(str(reaction_id)+'_obj', pathway_id)
        flux = float(flux_objective.getAnnotation().getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
        return str(reaction_id)+'_obj', float(flux)


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
        self.convertToCobra()
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
        self.convertToCobra()
        cobra_results = pfba(self.cobraModel, fraction_of_optimum)
        writeAnalysisResults(fbc_obj, cobra_results, pathway_id)


    ## Optimise for a target reaction (TODO) while fixing a source reaction to the fraction of its optimum
    #
    #
    def fractionReaction(self, source_reaction, target_reaction, fraction_of_source=0.75, pathway_id='rp_pathway'):
        #retreive the biomass objective and flux results and set as maxima
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        source_obj_id, source_flux = self.findCreateObjective(source_reaction, pathway_id)
        target_obj_id, target_flux = self.findCreateObjective(target_reaction, pathway_id)
        #set bounds biomass as a fraction
        old_upper_bound, old_lower_bound = self.setReactionConstraint(source_reaction, 
                                                                      source_flux*fraction_of_source, 
                                                                      source_flux*fraction_of_source)
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(target_obj_id), 
                'Setting active objective '+str(target_obj_id))
        self.convertToCobra()
        fbc_obj = fbc_plugin.getObjective(target_obj_id)
        self._checklibSBML(fbc_obj, 'Getting the following objective '+str(target_obj_id))
        cobra_results = self.cobraModel.optimize()
        self.writeAnalysisResults(fbc_obj, cobra_results)
        #reset the bounds to the original values for the target
        old_upper_bound, old_lower_bound = self.setReactionConstraint(source_reaction,
                                                                      old_upper_bound,
                                                                      old_lower_bound)
        return True


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
     

    ## Set a given reaction's upper and lower bounds
    #
    # Sets the upper and lower bounds of a reaction. Note that if the numerical values passed
    # are not recognised, new parameters are created for each of them
    #
    def setReactionConstraint(self, reaction_id, upper_bound, lower_bound):
        reaction = self.rpsbml.model.getReaction(reaction_id)
        if not reaction:
            self.logger.error('Cannot find the reaction: '+str(reaction_id))
            return False
        reac_fbc = reaction.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        ########## upper bound #############
        upper_param = None
        old_upper_param = self.rpsbml.model.getParameter(reac_fbc.getUpperFluxBound()).value
        for parameter in self.rpsbml.model.getListOfParameters():
            if parameter.getValue()==upper_bound:
                upper_param = parameter
        if not upper_param:
            upper_param = self.rpsbml.model.createParameter()
            self._checklibSBML(upper_param, 'creating target parameter')
            if upper_bound.is_integer():
                self._checklibSBML(upper_param.setId('B_'+str(upper_bound)),
                    'setting target parameter ID')
            else:
                self._checklibSBML(upper_param.setId('B_'+str(upper_bound).split('.')[0]+'_'+str(upper_bound).split('.')[1]), 
                    'setting target parameter ID')
            self._checklibSBML(upper_param.setSBOTerm(626),
                'setting target parameter SBO')
            #WARNING: need to have a way to generate the right unit
            self._checklibSBML(upper_param.setUnits('mmol_per_gDW_per_hr'),
                'setting target parameter Units')
            self._checklibSBML(upper_param.setValue(upper_bound),
                'setting target parameter Value')
            self._checklibSBML(upper_param.setConstant(True),
                'setting target parameter ID')
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_param.getId()), 
            'setting '+str(reaction_id)+' upper flux bound')
        ######### lower bound #############
        lower_param = None
        old_lower_param = self.rpsbml.model.getParameter(reac_fbc.getLowerFluxBound()).value
        for parameter in self.rpsbml.model.getListOfParameters():
            if parameter.getValue()==upper_bound:
                lower_param = parameter
        if not lower_param:
            lower_param = target_rpsbml.model.createParameter()
            self._checklibSBML(lower_param, 'creating target parameter')
            if upper_bound.is_integer():
                self._checklibSBML(lower_param.setId('B_'+str(upper_bound)), 
                    'setting target parameter ID')
            else:
                self._checklibSBML(lower_param.setId('B_'+str(upper_bound).split('.')[0]+'_'+str(upper_bound).split('.')[1]), 
                    'setting target parameter ID')
            self._checklibSBML(lower_param.setSBOTerm(626),
                'setting target parameter SBO')
            #WARNING: need to have a way to generate the right unit
            self._checklibSBML(lower_param.setUnits('mmol_per_gDW_per_hr'),
                'setting target parameter Units')
            self._checklibSBML(lower_param.setValue(upper_bound),
                'setting target parameter Value')
            self._checklibSBML(lower_param.setConstant(True),
                'setting target parameter ID')
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_param.getId()), 
            'setting '+str(reaction_id)+' lower flux bound')
        return old_upper_param, old_lower_param




                

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
