#!/usr/bin/env python3
"""
Add attributes and/or new child elements to XML files.

Basic usage examples:
  # see description in main
"""
from __future__ import annotations
import argparse
from logging import root
import os
from posixpath import split
import xml.etree.ElementTree as ET
import sys
import textwrap
from typing import List, Tuple
try:
    from lxml import etree as lxml_etree  # optional XSD validation
    from lxml import etree
except Exception:
    lxml_etree = None

import xml_formatter


def pretty_print_file(tree, outpath: str | None, indent: str, wrap_width: int, attr_per_line: bool,doo_pretty: bool):

    """
    Accept xml.etree.ElementTree.ElementTree or lxml Element/ElementTree.
    Convert xml.etree trees to lxml before serializing to get pretty_print support.
    """
    # If it's an xml.etree.ElementTree, convert to bytes then parse with lxml
    try:
        # detect xml.etree ElementTree/Element
        import xml.etree.ElementTree as ET_local
        is_etree = isinstance(tree, ET_local.ElementTree) or isinstance(tree, ET_local.Element)
    except Exception:
        ET_local = None
        is_etree = False

    if is_etree:
        # get bytes from xml.etree and parse with lxml
        xml_bytes = ET_local.tostring(tree.getroot() if isinstance(tree, ET_local.ElementTree) else tree, encoding="utf-8")
        lxml_root = etree.fromstring(xml_bytes)
        # ensure explicit end tags for empty elements
        for el in lxml_root.iter():
            if len(el) == 0 and el.text is None:
                el.text = ""
        xml_bytes_out = etree.tostring(lxml_root, pretty_print=True, xml_declaration=True, encoding='UTF-8')
        # write bytes
        with open(outpath, "wb") as f:
            f.write(xml_bytes_out)
    else:
        # assume it's already an lxml element/tree; handle both Element and ElementTree
        try:
            if hasattr(tree, "getroot"):
                lroot = tree.getroot()
            else:
                lroot = tree
            for el in lroot.iter():
                if len(el) == 0 and el.text is None:
                    el.text = ""
            xml_bytes_out = etree.tostring(lroot, pretty_print=True, xml_declaration=True, encoding='UTF-8')
            with open(outpath, "wb") as f:
                f.write(xml_bytes_out)
        except Exception:
            # fallback: use xml.etree write (ensures no <tag/> short form)
            try:
                if hasattr(tree, "write"):
                    tree.write(outpath, encoding='utf-8', xml_declaration=True, short_empty_elements=False)
                else:
                    # last resort: convert to string and write
                    s = ET_local.tostring(tree, encoding='utf-8')
                    with open(outpath, "wb") as f:
                        f.write(s)
            except Exception as e:
                raise
    xml_formatter.format_file(outpath)
    # 

def validate_with_xsd(tree: ET.ElementTree, xsd_path: str) -> Tuple[bool, str]:
    if not lxml_etree:
        return False, "lxml not installed"
    try:
        xml_bytes = ET.tostring(tree.getroot(), encoding="utf-8")
        xml_doc = lxml_etree.fromstring(xml_bytes)
        xsd_doc = lxml_etree.parse(xsd_path)
        schema = lxml_etree.XMLSchema(xsd_doc)
        ok = schema.validate(xml_doc)
        if ok:
            return True, ""
        else:
            return False, str(schema.error_log)
    except Exception as e:
        return False, str(e)

def update_constraint_element(xml_file,add_we,estimatorSolves):
    delete_old_schema=True
    with open(xml_file, 'rb') as xml_file:
     tree = ET.parse(xml_file)
    root = tree.getroot()
    compoFluidModel = True
    isThermal=False
    for elem in root.iter():
        #print(f"Tag: {elem.tag}, Attributes: {elem.attrib}, Text: {elem.text}")

        if elem.tag == "CompositionalMultiphaseWell" or elem.tag == "SinglePhaseWell":
            if elem.tag == "SinglePhaseWell":
                compoFluidModel = False
            #if "writeCSV" not in elem.attrib:
            #    esolves = '"'+str(estimatorSolves)+'"'
            #    elem.set("writeCSV","1")
            if add_we:
                 
                nlsTag="NonlinearSolverParameters"
                nlsAttributes={}
                nlsAttributes["newtonTol"]="1.0e-8"
                nlsAttributes["lineSearchAction"]="None"
                nlsAttributes["newtonMaxIter"]="20"
                elem.insert(0,ET.Element(nlsTag,nlsAttributes))
                lsTag="LinearSolverParameters"
                lsAttributes={}
                lsAttributes["directParallel"]="0"
                elem.insert(0,ET.Element(lsTag,lsAttributes))

            if  "isThermal" in elem.attrib:
                if elem.attrib["isThermal"] == "1":
                    isThermal = True
        if elem.tag == "WellControls":
            #print(f"Tag: {elem.tag}, Attributes: {elem.attrib}, Text: {elem.text}")
            #if "control" in elem.attrib:
            #    if delete_old_schema:
            #        elem.attrib.pop("control")
            if add_we:
                elem.set("estimateWellSolution",str(estimatorSolves))
            isProducer = elem.attrib['type'] == 'producer'
            if isProducer:
                constraintType="Production"
                pressureType="Minimum"
            else:
                constraintType="Injection"
                pressureType="Maximum"
            # setup phase constraint
            if 'targetPhaseName' in elem.attrib:
                writePhaseConstraint=True
                phaseConstraintTag= constraintType+"PhaseVolumeRateConstraint"
                phaseConstraintAttributes={}
                phaseConstraintAttributes["name"]="max"+elem.attrib['targetPhaseName'].lower() +"prod" if isProducer else "inj"
                if 'targetPhaseName' in elem.attrib:
                    phaseConstraintAttributes["phaseName"] = elem.attrib['targetPhaseName']
                    if delete_old_schema:
                        elem.attrib.pop("targetPhaseName")
                if "targetPhaseRateTableName" in elem.attrib:
                    phaseConstraintAttributes["constraintScheduleTableName"]=elem.attrib["targetPhaseRateTableName"]
                    if delete_old_schema:
                        elem.attrib.pop("targetPhaseRateTableName")
                elif "targetPhaseRate" in elem.attrib:
                    phaseConstraintAttributes["phaseRate"]=elem.attrib["targetPhaseRate"]
                    if delete_old_schema:
                        elem.attrib.pop("targetPhaseRate")
                else:
                    writePhaseConstraint=False
                    print("error missing phase rate info")
                if not isProducer:
                    if compoFluidModel:
                        if "injectionStream" in elem:
                            phaseConstraintAttributes["injectionStream"]=elem.attrib["injectionStream"]
                            phaseConstraintAttributes["injectionTemperature"]=elem.attrib["injectionTemperature"]
                            if delete_old_schema:
                                elem.attrib.pop("injectionStream")
                                elem.attrib.pop("injectionTemperature") 
                        else:
                            print("error missinging injectionStream ",elem.attrib)

                if writePhaseConstraint:
                    elem.append(ET.Element(phaseConstraintTag,phaseConstraintAttributes))

            # setup pressure constraint
            pressureConstraintTag= pressureType+"BHPConstraint"
            pressureConstraintAttributes={}
            pressureConstraintAttributes["name"]="minbhp" if isProducer else "maxbhp"
            if 'targetBHPTableName' in elem.attrib:
                pressureConstraintAttributes["constraintScheduleTableName"]=elem.attrib["targetBHPTableName"]
                if delete_old_schema:
                    elem.attrib.pop("targetBHPTableName")   
            elif 'targetBHP' in elem.attrib:
                pressureConstraintAttributes["targetBHP"]=elem.attrib["targetBHP"]
                if delete_old_schema:
                    elem.attrib.pop("targetBHP")    
            else:
                print('error missing bhp info')

            if 'referenceElevation' in elem.attrib:
                pressureConstraintAttributes["referenceElevation"]=elem.attrib["referenceElevation"]
                if delete_old_schema:
                    elem.attrib.pop("referenceElevation")   
            else:
                print('error missing bhp referenceElevation')

            elem.append(ET.Element(pressureConstraintTag,pressureConstraintAttributes))
 
            totalVolRateTag= constraintType+"VolumeRateConstraint"
            totalVolRateAttributes={}
            totalVolRateAttributes["name"]="maxvolrateprod" if isProducer else  "maxvolrateinj"
            if 'targetTotalRate' in elem.attrib or  'targetTotalRateTableName' in elem.attrib:
                if 'targetTotalRate' in elem.attrib:
                    totalVolRateAttributes["volumeRate"]=elem.attrib["targetTotalRate"]
                    if delete_old_schema:
                        elem.attrib.pop("targetTotalRate")  
                if 'targetTotalRateTableName' in elem.attrib:
                    totalVolRateAttributes["constraintScheduleTableName"]=elem.attrib["targetTotalRateTableName"]
                    if delete_old_schema:
                        elem.attrib.pop("targetTotalRateTableName")
                if not isProducer:
                    if compoFluidModel:
                        totalVolRateAttributes["injectionStream"]=elem.attrib["injectionStream"]
                        totalVolRateAttributes["injectionTemperature"]=elem.attrib["injectionTemperature"]
                        if delete_old_schema:
                            elem.attrib.pop("injectionStream")
                            elem.attrib.pop("injectionTemperature")
                    else:
                        if isThermal:
                            if "injectionStream" in elem.attrib:
                                totalVolRateAttributes["injectionStream"]=elem.attrib["injectionStream"]
                                totalVolRateAttributes["injectionTemperature"]=elem.attrib["injectionTemperature"]
                                if delete_old_schema:
                                    elem.attrib.pop("injectionStream")
                                    elem.attrib.pop("injectionTemperature")
                            else:
                                totalVolRateAttributes["injectionStream"]="{1.0}"
                                totalVolRateAttributes["injectionTemperature"]="123" 
                
                elem.append(ET.Element(totalVolRateTag,totalVolRateAttributes))

            totalMassRateTag= constraintType+"MassRateConstraint"
            totalMassRateAttributes={}
            totalMassRateAttributes["name"]="maxmassrateprod" if isProducer else  "maxmassrateinj"
            if 'targetMassRate' in elem.attrib or  'targetMassRateTableName' in elem.attrib:
                if 'targetMassRate' in elem.attrib:
                    totalMassRateAttributes["massRate"]=elem.attrib["targetMassRate"]
                    if delete_old_schema:
                        elem.attrib.pop("targetMassRate")
                if 'targetMassRateTableName' in elem.attrib:
                    totalMassRateAttributes["constraintScheduleTableName"]=elem.attrib["targetMassRateTableName"]
                    if delete_old_schema:
                        elem.attrib.pop("targetMassRateTableName")
                if not isProducer:
                    if True: #compoFluidModel:
                        totalMassRateAttributes["injectionStream"]=elem.attrib["injectionStream"]
                        totalMassRateAttributes["injectionTemperature"]=elem.attrib["injectionTemperature"]
                        if delete_old_schema:
                            elem.attrib.pop("injectionStream")
                            elem.attrib.pop("injectionTemperature")
                elem.append(ET.Element(totalMassRateTag,totalMassRateAttributes))

    wc={}
    isCompositional=False
    nlsparms={}
    useMass=False
    useMasVal="0"
    for elem in root.iter():
        if elem.tag == "CompositionalMultiphaseWell": 
            isCompositional=True
            if "maxRelativePressureChange" in elem.attrib:
                nlsparms["maxRelativePressureChange"]=elem.attrib["maxRelativePressureChange"]
                elem.attrib.pop("maxRelativePressureChange")
            if "maxCompFractionChange" in elem.attrib:
                nlsparms["maxCompFractionChange"]=elem.attrib["maxCompFractionChange"]
                elem.attrib.pop("maxCompFractionChange")
            if "useMass" in elem.attrib:
                useMass=True
                useMassVal=elem.attrib["useMass"]
                elem.attrib.pop("useMass")
            if "writeCSV" in elem.attrib:
                nlsparms["writeCSV"]=elem.attrib["writeCSV"]
                elem.attrib.pop("writeCSV")    
 
            break
        elif elem.tag == "SinglePhaseWell":
            if "writeCSV" in elem.attrib:
                nlsparms["writeCSV"]=elem.attrib["writeCSV"]
                elem.attrib.pop("writeCSV") 
    for elem in root.iter():
        if elem.tag == "CompositionalMultiphaseWell" or elem.tag == "SinglePhaseWell":
            elem.tag="WellManager"
            if useMass:
                elem.set("useMass",useMassVal)
        
        if elem.tag == "WellControls":
            for k in nlsparms:
                elem.attrib[k] = nlsparms[k]
            if isCompositional:
                elem.tag="CompositionalMultiphaseWell"
            else:    
                elem.tag="SinglePhaseWell"

    return True ,tree  

def removeUseMass(xml_file,delete_old_schema,add_we,estimatorSolves):
    with open(xml_file, 'rb') as xml_file:
     tree = ET.parse(xml_file)
    root = tree.getroot()
    compoFluidModel = True
    isThermal=False

    wc={}
    isCompositional=False
    nlsparms={}
    useMass=False
    useMasVal="0"
    mod=False
    for elem in root.iter():
        if elem.tag == "CompositionalMultiphaseWell": 
            isCompositional=True
            if "useMass" in elem.attrib:
                elem.attrib.pop("useMass")
                mod=True

    return mod,tree  

def add_we(xml_file,estimatorSolves):
    with open(xml_file, 'rb') as xml_file:
     tree = ET.parse(xml_file)
    root = tree.getroot()
    compoFluidModel = True
    isThermal=False

    wc={}
    isCompositional=False
    nlsparms={}
    useMass=False
    useMasVal="0"
    mod=True
    for elem in root.iter():
        if elem.tag == "CompositionalMultiphaseWell" or elem.tag == "SinglePhaseWell": 
            elem.attrib["estimateWellSolution"] = str(estimatorSolves)

    return mod,tree  

def main1(ifs,ofs,add_we,estimatorSolves,args):    
    try:
        mod,tree= update_constraint_element(ifs, add_we,estimatorSolves)
        #mod,tree= removeUseMass(ifs, delete_old_schema,add_we,estimatorSolves)
        if mod:
            wrap_width=80
            indent="  "
            attr_per_line=True
            pretty_print_file(tree, ofs,  indent,  wrap_width,  attr_per_line,False)
            # Optionally validate against XSD
            if args.xsd:
                ok, msg = validate_with_xsd(tree, args.xsd)
                if not ok:
                    print("XSD validation failed:", msg, file=sys.stderr)
                    # still write a diagnostic copy
                    tree.write(ofs + ".failed.xml", encoding="utf-8", xml_declaration=True)
                    print("XSD failed for ",ofs )
                    xml_formatter.format_file(ofs + ".failed.xml")
                else:
                    print("XSD validation passed. Wrote to ",ofs)
    except Exception as e:
        print("Error occurred:", e,ifs)
        
def add_estimator(ifs,ofs,estimatorSolves,args):    
    try:
        #tree= updateconstraint_element(ifs, delete_old_schema,add_we,estimatorSolves)
        mod,tree= add_we(ifs, estimatorSolves)
        if mod:
            wrap_width=80
            indent="  "
            attr_per_line=True
            pretty_print_file(tree, ofs,  indent,  wrap_width,  attr_per_line,args.pretty)
            # Optionally validate against XSD
            if args.xsd:
                ok, msg = validate_with_xsd(tree, args.xsd)
                if not ok:
                    print("XSD validation failed:", msg, file=sys.stderr)
                    # still write a diagnostic copy
                    tree.write(ofs + ".failed.xml", encoding="utf-8", xml_declaration=True)
                    print("XSD failed for ",ofs )
                    xml_formatter.format_file(ofs + ".failed.xml")
                else:
                    print("XSD validation passed. Wrote to ",ofs)
    except Exception as e:
        print("Error occurred:", e,ifs)

def main():
    # Two modes: single file or batch from list
    # To create a list of files to process
    #   egrep -l -m1 -R "WellControl|SinglePhaseWell|CompositionalMultiphaseWell" ./ --include=\*.xml > convert_lst.txt
    ## convert from file list 
    # cws_v3 -d --xsd=/Users/byer3/GEOS-DEV-1105/wc1014/src/coreComponents/schema/schema.xsd -r -f convert_lst.txt | tee convert.log
    ## convert single file inplace
    # cws_v3.py -scompositionalMultiphaseWell/compositional_multiphase_wells_1d.xml --xsd=/Users/byer3/GEOS-DEV/wmwe0603/src/coreComponents/schema/schema.xsd 
 
    p = argparse.ArgumentParser(description="Add attributes and/or child elements to XML")
    p.add_argument("-s", "--sourcefile", required=False, help="source XML file")
    p.add_argument("-t", "--targetfile", required=False, help="output XML file")
    p.add_argument("-a", "--add",  action="store_true",default=False,  help="add we strings")
    p.add_argument("-e", "--estimatorSolves",  type=int ,default=0, help="when to use estimator")
    p.add_argument("-r", "--replace", action="store_true",default=True,    help="in place substitution")
    p.add_argument("-f", "--file",  type=str ,default="", help="file with list of files to process")
    p.add_argument("--xsd", required=True, help=" XSD file to validate output (requires lxml)")
    args = p.parse_args()

    ofn_tag=""
    if args.estimatorSolves >0 and not args.add:
        ofn_tag+="we"+str(args.estimatorSolves)
    if args.file:
        if os.path.exists(args.file):
            ifs = open(args.file,"r")
            for f in ifs:           
                args.sourcefile=f.rstrip().lstrip()
                if args.replace:
                    fn1=f.rstrip().lstrip()
                    fn2=fn1
                else:
                    fns = f.split()
                    fn1=fns[0].rstrip().lstrip()
                    if len(ofn_tag)>0:
                        fn2=fns[0].rstrip().lstrip()
                        fn2=fn2.replace(".xml",ofn_tag+".xml")
                    else:
                        fn2=fns[1].rstrip().lstrip()
                print("Processing ",fn1,fn2,args.add,args.estimatorSolves)
                if args.estimatorSolves >0 :
                    add_estimator(fn1,fn2,args.estimatorSolves,args)
                else:
                    main1(fn1,fn2,args.add,args.estimatorSolves,args)
        else:
            print("File with list of files to process not found ",args.file)
        
    else:
        if os.path.exists(args.sourcefile):
             if not args.targetfile:
                args.targetfile=args.sourcefile
             main1(args.sourcefile,args.targetfile,args.add,args.estimatorSolves,args)
        else:
            print("Source file not found ",args.sourcefile)

if __name__ == "__main__":
    main()