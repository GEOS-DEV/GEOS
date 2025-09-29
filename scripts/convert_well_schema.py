import xml.etree.ElementTree as ET
import os
import optparse 

def print_xml_with_indent_and_attributes(ofn,tree, indent="  "):
    import textwrap
    ofs = open(ofn,"w")
    root = tree.getroot()
    ET.indent(tree, space=indent, level=0)  # Indent nested elements

    def format_element(elem, level=0):
        spaces = indent * level
        tag_line = f"{spaces}<{elem.tag}"
        attrib_lines = []
        for k, v in elem.attrib.items():
            v =v.replace(" ","")
            v = v.replace(",",", ")
            attr_str = f'{k}="{v}"'
            wrapped = textwrap.wrap(attr_str, width=100, break_long_words=False, break_on_hyphens=False)
            for i, line in enumerate(wrapped):
                if i == 0:
                    attrib_lines.append(f"{indent * (level + 1)}{indent}{line}")
                else:
                    # Extra indent for all rows after the first
                    attrib_lines.append(f"{indent * (level + 2)}{"  "+line}")
        if attrib_lines:
            tag_line += "\n" + "\n".join(attrib_lines)
            tag_line += f">"
        else:
            istrt = tag_line.find("<") 
            if istrt > -1:
                tag_line += ">"
            else:
                tag_line += ">"
        lines = [tag_line]
        if elem.text and elem.text.strip():
            lines.append(f"{spaces}{indent}{elem.text.strip()}")
        for child in elem:
            lines.append(format_element(child, level + 1))
        lines.append(f"{spaces}</{elem.tag}>")
        return "\n".join(lines)
    ofs.write(format_element(root))
    ofs.close()

def update_or_add_constraint_element(xml_file, ofs):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    for elem in root.iter():
        #print(f"Tag: {elem.tag}, Attributes: {elem.attrib}, Text: {elem.text}")
        if elem.tag == "CompositionalMultiphaseWell":
            if "writeCSV" not in elem.attrib:
                elem.set("writeCSV","1")
            elem.set("useNewCode", "1")
            nlsTag="NonlinearSolverParameters"
            nlsAttributes={}
            nlsAttributes["newtonTol"]="1.0e-8"
            nlsAttributes["lineSearchAction"]="None"
            nlsAttributes["newtonMaxIter"]="20"
            elem.insert(0,ET.Element(nlsTag,nlsAttributes))
            lsTag="LinearSolverParamters"
            lsAttributes={}
            lsAttributes["directParallel"]="0"
            elem.insert(0,ET.Element(lsTag,lsAttributes))
            
        if elem.tag == "WellControls":
            #print(f"Tag: {elem.tag}, Attributes: {elem.attrib}, Text: {elem.text}")
            elem.set("estimateWellSolution","2")
            isProducer = elem.attrib['type'] == 'producer'
            if isProducer:
                constraintType="Production"
                pressureType="Minimum"
            else:
                constraintType="Injection"
                pressureType="Maximum"
            # setup phase constraint
            if 'targetPhaseName' in elem.attrib:
                phaseConstraintTag= "Phase"+constraintType+"Constraint"
                phaseConstraintAttributes={}
                phaseConstraintAttributes["name"]="max"+elem.attrib['targetPhaseName'].lower() +"prod" if isProducer else "max"+elem.attrib['targetPhaseName'].lower() +"inj"
                if 'targetPhaseName' in elem.attrib:
                    phaseConstraintAttributes["phaseName"] = elem.attrib['targetPhaseName']
                if "targetPhaseRateTableName" in elem.attrib:
                    phaseConstraintAttributes["constraintScheduleTableName"]=elem.attrib["targetPhaseRateTableName"]
                elif "targetPhaseRate" in elem.attrib:
                    phaseConstraintAttributes["phaseRate"]=elem.attrib["targetPhaseRate"]
                else:
                    print("error missing phase rate info")
                if not isProducer:
                    phaseConstraintAttributes["injectionStream"]=elem.attrib["injectionStream"]
                    phaseConstraintAttributes["injectionTemperature"]=elem.attrib["injectionTemperature"]
                
                elem.append(ET.Element(phaseConstraintTag,phaseConstraintAttributes))

            # setup pressure constraint
            pressureConstraintTag= pressureType+"BHPConstraint"
            pressureConstraintAttributes={}
            pressureConstraintAttributes["name"]="minbhp" if isProducer else "maxbhp"
            if 'targetBHPTableName' in elem.attrib:
                pressureConstraintAttributes["constraintScheduleTableName"]=elem.attrib["targetBHPTableName"]
            elif 'targetBHP' in elem.attrib:
                pressureConstraintAttributes["targetBHP"]=elem.attrib["targetBHP"]
            else:
                print('error missing bhp info')

            if 'referenceElevation' in elem.attrib:
                pressureConstraintAttributes["referenceElevation"]=elem.attrib["referenceElevation"]
            else:
                print('error missing bhp referenceElevation')

            elem.append(ET.Element(pressureConstraintTag,pressureConstraintAttributes))
 
            totalVolRateTag= "TotalVol"+constraintType+"Constraint"
            totalVolRateAttributes={}
            if 'targetTotalRate' in elem.attrib or  'targetTotalRateTableName' in elem.attrib:
                if 'targetTotalRate' in elem.attrib:
                    totalVolRateAttributes["volumeRate"]=elem.attrib["targetTotalRate"]
                if 'targetTotalRateTableName' in elem.attrib:
                    totalVolRateAttributes["constraintScheduleTableName"]=elem.attrib["targetTotalRateTableName"]
                if not isProducer:
                    totalVolRateAttributes["injectionStream"]=elem.attrib["injectionStream"]
                    totalVolRateAttributes["injectionTemperature"]=elem.attrib["injectionTemperature"]
                elem.append(ET.Element(totalVolRateTag,totalVolRateAttributes))
        print_xml_with_indent_and_attributes(ofs,tree)

 
def main1(ifs,ofs):    
    update_or_add_constraint_element(ifs,ofs)


if __name__ == "__main__":
 
    bdir = "/Users/byer3/geos_models/geos-total-dataset/GreatNorthernLight/GNL_FlowOnly/depletion"
    #fn = "/Users/byer3/opm/opm-simulators/tests/include/b1_vfp_flowline.inc"
    ifn = 'GNL_BO_WELL_FULLRATE.xml'
    ofn = 'GNL_BO_WELL_FULLRATE_WE2.xml'
    #ofn="/Users/byer3/GEOS-DEV-1105/we0708/inputFiles/compositionalMultiphaseWell/include"
    #main1(os.path.join(bdir,ifn),os.path.join(bdir,ofn))

    parser = optparse.OptionParser()
    parser.add_option("-s", "--sourcefile",  type="str" ,default="", help="source file")
    parser.add_option("-t", "--targetfile",  type="str" ,default="", help="target file")
    (options, args) = parser.parse_args()
    main1(options.sourcefile,options.targetfile)



    print('Finished')