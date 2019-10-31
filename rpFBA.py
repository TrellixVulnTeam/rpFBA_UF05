import cobra
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
            #add the results of FBA run to the annotation of FBA objective
            for flux_obj in fbc_obj.getListOfFluxObjectives():
                obj_annot = flux_obj.getAnnotation()
                brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:obj units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(flux_obj.getReaction()))+'" /> </brsynth:brsynth>')
                brsynth_annot.addChild(tmpAnnot.getChild('obj'))
            #add the results of the FBA reactions for each rp_pathway reactions 
            #TODO: need to take care of the case where the annotation exists already --> overwrite
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                reac_annot = reac.getAnnotation()
                brsynth_annot = reac_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:fba_'+str(fbc_obj.getId())+' units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(member.getIdRef()))+'" /> </brsynth:brsynth>')
                brsynth_annot.addChild(tmpAnnot.getChild('fba_'+str(fbc_obj.getId())))
            del res
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


    ## function to add the splitObj or any other objective that does not exist in the model
    #
    #
    def addObjective():
        pass       


    ## Given that there can be multiple objectives defined in the SBML, this function switches between different ones
    #
    # TODO: 
    def switchObjective(self, objId):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        listObj = fbc_plugin.getListOfObjectives()
        if objId in [i.getId() for i in listObj]:
            listObj.setActiveObjective(objId)
            return objId
        else:
            logger.warning('The objective Id '+str(objId)+' does not exist in rpsbml')
            return None
        self.libsbml_to_cobra()


    ##
    #
    #
    def libsbml_to_cobra(self):
        self.cobraModel = cobra.io.read_sbml_model(document.toXMLNode().toXMLString())
        #for an old version of cobrapy (0.4)
        #return cobra.io.sbml3.parse_xml_into_model(lxml.etree.fromstring(rpsbml.document.toXMLNode().toXMLString()))
        
    
    ##
    #
    #TODO: should update the 
    def simulate(self):
        res = self.cobraModel.optimize()
        return self.cobraModel.summary()
        
    
    ##
    #
    # child: rule_score, fba_biomass_score, fba_target_score, fba_splitObj_score, global_score
    def updateFBAscore(self, child, value):
        self.rpsbml()
        #reaction
        for reac in self.rpsbml.getListOfReactions():
            reac.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(child).addAttr('value', str(value))
        #pathway


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

