#!/usr/bin/env python3
"""Extract source/schema metadata used by the GEOS-MPM manual.

The output JSON is an intermediate generated artifact consumed by
tools/render_generated_tex.py.  This script does not write any human-maintained
LaTeX source.
"""
import argparse, re, json, pathlib, ast, hashlib
from collections import defaultdict, Counter

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=pathlib.Path, required=True, help="path to the GEOS checkout root")
parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path("generated/data/geos_mpm_extracted.json"), help="metadata JSON output path")
parser.add_argument("--write-schema-docs", action="store_true", help="refresh MPM wrapper reference RST files")
args = parser.parse_args()
repo = args.repo.resolve()

def read(p):
    return (repo/p).read_text(errors='ignore')

def line_no(text, idx):
    return text.count('\n',0,idx)+1

def clean_cpp_string_arg(arg):
    # concatenate C++ string literals, drop computed expressions and enum concat, convert escaped whitespace
    parts = re.findall(r'"((?:\\.|[^"\\])*)"', arg, flags=re.S)
    if not parts:
        return ''
    s=''.join(parts)
    s=s.replace('\\n','; ').replace('\\t',' ')
    s=s.replace('\\"','"')
    s=re.sub(r'\s+', ' ', s).strip()
    # clean punctuation from enum concatenation leftovers
    return s

def view_key_map(text):
    m={}
    for name,val in re.findall(r'static\s+constexpr\s+char\s+const\s*\*\s*([A-Za-z0-9_]+)String\s*\(\s*\)\s*\{\s*return\s+"([^"]+)"\s*;\s*\}', text):
        m[f'viewKeyStruct::{name}String()']=val
    return m

def extract_register_wrappers(path, header_paths=[]):
    text=read(path)
    # Resolve simple constructor defaults without importing or running GEOS.
    ctor_text=text[:text.find('registerWrapper')]
    ctor_defaults=dict(re.findall(r'\b(m_\w+)\s*\(\s*([^()]*)\s*\)\s*[,\n]', ctor_text))
    # Recover nested initializer expressions such as numeric_limits::max().
    for match in re.finditer(r'\b(m_\w+)\s*\(', ctor_text):
        start=match.end(); end=start; depth=1
        while end < len(ctor_text) and depth:
            if ctor_text[end]=='(':depth+=1
            elif ctor_text[end]==')':depth-=1
            end+=1
        if depth==0:ctor_defaults[match[1]]=ctor_text[start:end-1].strip()
    member_types={}
    for header in header_paths:
        member_types.update({name:typ.strip() for typ,name in re.findall(
            r'^\s*([\w:]+(?:\s*<[^;]+?>)?)\s+(m_\w+)\s*;', read(header), re.M)})
    maps={}
    for h in header_paths:
        maps.update(view_key_map(read(h)))
    maps.update(view_key_map(text))
    # manual scanner: find registerWrapper(, then collect statement to ; at top level wrt quotes/parens roughly
    rows=[]
    i=0
    while True:
        idx=text.find('registerWrapper', i)
        if idx<0: break
        # ensure not registerWrapperBase? okay if registerWrapper<...> included maybe skip template wrapper on node fields later by path
        j=idx
        in_str=False; esc=False; par=0; started=False
        while j<len(text):
            c=text[j]
            if in_str:
                if esc: esc=False
                elif c=='\\': esc=True
                elif c=='"': in_str=False
            else:
                if c=='"': in_str=True
                elif c=='(': par+=1; started=True
                elif c==')': par=max(0,par-1)
                elif c==';' and started and par==0:
                    j+=1; break
            j+=1
        stmt=text[idx:j]
        i=j
        # Skip nodeManager.registerWrapper<...> and nodeSets.registerWrapper in solver field registrations; retain only direct solver registerWrapper calls.
        prefix=text[max(0, idx-40):idx]
        if '.' in prefix.split('\n')[-1] or '->' in prefix.split('\n')[-1]:
            continue
        m=re.search(r'registerWrapper\s*\(\s*([^,]+),', stmt, flags=re.S)
        if not m: continue
        raw=m.group(1).strip()
        if raw.startswith('"'):
            name=re.match(r'"([^"]+)"', raw).group(1)
        else:
            name=maps.get(raw, raw)
        # input flag
        flag=''
        mf=re.search(r'setInputFlag\s*\(\s*InputFlags::([A-Z_]+)\s*\)', stmt)
        if mf: flag=mf.group(1)
        # restart
        restart=''
        mr=re.search(r'setRestartFlags\s*\(\s*RestartFlags::([A-Z_]+)\s*\)', stmt)
        if mr: restart=mr.group(1)
        default=''
        md=re.search(r'set(?:Apply)?DefaultValue\s*\((.*?)\)', stmt, flags=re.S)
        if md:
            default=re.sub(r'\s+', ' ', md.group(1).strip())
        # descriptions can have nested ) from EnumStrings; capture between setDescription( and following ); approximated
        desc=''
        pos=stmt.find('setDescription')
        if pos>=0:
            # find opening paren
            op=stmt.find('(', pos)
            k=op+1; par2=1; in_s=False; esc2=False
            while k<len(stmt):
                c=stmt[k]
                if in_s:
                    if esc2: esc2=False
                    elif c=='\\': esc2=True
                    elif c=='"': in_s=False
                else:
                    if c=='"': in_s=True
                    elif c=='(': par2+=1
                    elif c==')':
                        par2-=1
                        if par2==0: break
                k+=1
            desc=clean_cpp_string_arg(stmt[op+1:k])
        default_expression=default
        if default in ctor_defaults:
            default=ctor_defaults[default].strip() or '{}'
        enum_match=re.fullmatch(r'mpm::(\w+)::(\w+)', default)
        if enum_match:
            enum_text=read(pathlib.Path('src/coreComponents/physicsSolvers/solidMechanics/MPMSolverEnums.hpp'))
            body=re.search(r'enum\s+(?:class|struct)\s+'+enum_match[1]+r'[^\{]*\{(.*?)\};',enum_text,re.S)
            strings=re.search(r'ENUM_STRINGS\(\s*mpm::'+enum_match[1]+r'\s*,(.*?)\);',enum_text,re.S)
            if body and strings:
                symbols=[v.strip().split('=')[0].strip() for v in re.sub(r'//[^\n]*','',body[1]).split(',') if v.strip()]
                labels=re.findall(r'"([^"]+)"',strings[1])
                if enum_match[2] in symbols and len(labels)==len(symbols):default=labels[symbols.index(enum_match[2])]
        member=re.search(r'registerWrapper\s*\([^,]+,\s*&\s*(m_\w+)',stmt,re.S)
        typ=member_types.get(member[1], '') if member else ''
        if typ.startswith('mpm::'):typ='string (enum)'
        rows.append({'name':name,'flag':flag,'default':default,'default_expression':default_expression,
                     'type':typ,'restart':restart,'description':desc,'source':str(path),'line':line_no(text,idx)})
    return rows

# Solver wrappers
solver_wrappers=extract_register_wrappers(pathlib.Path('src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp'), [pathlib.Path('src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp')])
# Keep only direct solver constructor wrappers before registerDataOnMesh line; scanner already skips member nodeManager. But include all direct.
# Events wrappers
base_wrappers=extract_register_wrappers(pathlib.Path('src/coreComponents/events/mpmEvents/MPMEventBase.cpp'), [pathlib.Path('src/coreComponents/events/mpmEvents/MPMEventBase.hpp')])
event_data=[]
for h in sorted((repo/'src/coreComponents/events/mpmEvents').glob('*MPMEvent.hpp')):
    if h.name in ('MPMEventBase.hpp','MPMEvents.hpp'): continue
    cpp=h.with_suffix('.cpp')
    if not cpp.exists(): continue
    htext=h.read_text(errors='ignore')
    cat=re.search(r'static\s+string\s+catalogName\s*\(\s*\)\s*\{\s*return\s+"([^"]+)"', htext)
    cname=cat.group(1) if cat else h.stem.replace('MPMEvent','')
    doc=''
    mdoc=re.search(r'\* @class.*?\n \*\s*\n \*\s*(.*?)\n \*/', htext, flags=re.S)
    if mdoc:
        doc=re.sub(r'\n\s*\*\s*',' ',mdoc.group(1)).strip()
    wrappers=base_wrappers + extract_register_wrappers(pathlib.Path(cpp.relative_to(repo)), [pathlib.Path(h.relative_to(repo)), pathlib.Path('src/coreComponents/events/mpmEvents/MPMEventBase.hpp')])
    event_data.append({'name':cname,'class':h.stem,'description':doc,'attributes':wrappers,'source':str(h.relative_to(repo))})

# Particle types and columns
def extract_enum_strings(path, enum_qual):
    text=read(path)
    pat=re.escape(enum_qual)+r'\s*,(.*?)\);'
    m=re.search(pat, text, flags=re.S)
    if not m: return []
    return re.findall(r'"([^"]+)"', m.group(1))
particle_types=extract_enum_strings(pathlib.Path('src/coreComponents/mesh/ParticleType.hpp'), 'ParticleType')
particle_columns=extract_enum_strings(pathlib.Path('src/coreComponents/mesh/particleGenerators/ParticleMeshGenerator.hpp'), 'ParticleMeshGenerator::ParticleColumnHeaders')
# Particle header names, not positional indices, select columns. Parse the
# actual default switch, including fields previously omitted by the hardcoded list.
reader_text=read(pathlib.Path('src/coreComponents/mesh/particleGenerators/ParticleMeshGenerator.cpp'))
switch=reader_text[reader_text.index('switch( static_cast< ParticleColumnHeaders >( c ) )'):]
switch=switch[:switch.index('default:')]
col_defaults={}
for group in switch.split('break;'):
    value=re.search(r'defaultValue\s*=\s*([^;]+);',group)
    if value:
        for name in re.findall(r'case ParticleColumnHeaders::(\w+):',group):col_defaults[name]=value[1].strip()
particle_columns_data=[]
for index,name in enumerate(particle_columns):
    geometry=name in ('GeometryVersion','GeometryShape','ParticleVolume') or name.startswith('Geometry')
    particle_columns_data.append({'index':index,'name':name,
        'default':'required schema 2' if geometry else col_defaults.get(name,'no default'),
        'required':geometry or name not in col_defaults,
        'required_when':'3D only' if name in ('PositionZ','GeometryXZ','GeometryYZ','GeometryZX','GeometryZY','GeometryZZ') else ''})

# Solver enum options live primarily in MPMSolverEnums.hpp. Preserve serialized
# XML spellings, which need not match C++ enum names.
htext='\n'.join(read(pathlib.Path('src/coreComponents/physicsSolvers/solidMechanics')/name)
                for name in ['MPMSolverEnums.hpp','SolidMechanicsMPM.hpp'])
enums=[]
for m in re.finditer(r'enum\s+(?:class|struct)\s+([A-Za-z0-9_]+)[^{]*\{(.*?)\};', htext, flags=re.S):
    name=m.group(1)
    body=m.group(2)
    vals=[]
    for line in body.splitlines():
        line=line.split('//')[0].split('//!<')[0].strip().rstrip(',')
        if not line: continue
        line=line.split('=')[0].strip()
        if re.match(r'[A-Za-z_]', line): vals.append(line)
    strings=re.search(r'ENUM_STRINGS\(\s*mpm::'+name+r'\s*,(.*?)\);',htext,re.S)
    labels=re.findall(r'"([^"]+)"',strings[1]) if strings else []
    enums.append({'name':name,'values':labels or vals,'cpp_values':vals})

particle_region_wrappers=extract_register_wrappers(pathlib.Path('src/coreComponents/mesh/ParticleRegion.cpp'),
    [pathlib.Path('src/coreComponents/mesh/ParticleRegion.hpp')])

# PFW parameters table: extract dict literal region approximately and ast parse
pfw_text=read(pathlib.Path('scripts/preProcessing/materialPointMethodParticleFileWriter/particleFileWriter.py'))
params=[]
start=pfw_text.find('parameters = {')
if start>=0:
    # collect balanced braces ignoring strings/comments roughly: easiest use ast parse of assignment after replacing comments okay
    i=start+pfw_text[start:].find('{')
    j=i; brace=0; in_s=False; quote=''; esc=False
    while j<len(pfw_text):
        c=pfw_text[j]
        if in_s:
            if esc: esc=False
            elif c=='\\': esc=True
            elif c==quote: in_s=False
        else:
            if c in ('"',"'"): in_s=True; quote=c
            elif c=='{': brace+=1
            elif c=='}':
                brace-=1
                if brace==0:
                    j+=1; break
        j+=1
    literal=pfw_text[i:j]
    try:
        # Remove comments not necessary for ast; comments okay. Use literal_eval after wrapping? dict values can have None/False okay.
        parsed=ast.literal_eval(literal)
        for k,v in parsed.items():
            default=v[0]
            emit=v[1]
            params.append({'name':str(k),'default':repr(default),'emitsSolverAttribute':bool(emit)})
    except Exception as e:
        # fallback regex
        for m in re.finditer(r"['\"]([^'\"]+)['\"]\s*:\s*\((.*?),\s*(True|False)\s*\)", literal):
            params.append({'name':m.group(1),'default':m.group(2).strip(),'emitsSolverAttribute':m.group(3)=='True'})

# PFW grid field presets parse
silo_common=re.findall(r'_PFW_SILO_COMMON_GRID_FIELDS\s*=\s*\[(.*?)\]', pfw_text, re.S)
silo_all=re.findall(r'_PFW_SILO_ALL_GRID_FIELDS\s*=\s*\[(.*?)\]', pfw_text, re.S)
def list_strings(s): return re.findall(r'"([^"]+)"', s)
silo_common=list_strings(silo_common[0]) if silo_common else []
silo_all=list_strings(silo_all[0]) if silo_all else []

# PFW geometry classes
geom_text=read(pathlib.Path('scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_geometryObjects.py'))
geom_classes=[]
for m in re.finditer(r'^class\s+([A-Za-z_][A-Za-z0-9_]*)\s*(?:\(([^)]*)\))?:', geom_text, re.M):
    name=m.group(1); bases=m.group(2) or ''
    lineno=line_no(geom_text,m.start())
    # get docstring/description? many no docstrings. Categorize based on base/name.
    geom_classes.append({'name':name,'bases':bases,'line':lineno})
# separate object vs wrappers/utilities
geometry_names=[]; wrapper_names=[]; utility_classes=[]
for row in geom_classes:
    name=row['name']; bases=row['bases']
    if name in ('Geometry','BaseWrapper','SetOperation','Grid','maxFlawRadius','maxFlawRadius2D','MaterialProperties'):
        utility_classes.append(row)
    elif 'Wrapper' in name or name.endswith('Wrapper') or 'BaseWrapper' in bases or name in ('transform',):
        wrapper_names.append(row)
    elif 'Geometry' in bases or bases.strip()=='Geometry' or name in ('union','intersection','difference') or 'SetOperation' in bases:
        geometry_names.append(row)
    else:
        utility_classes.append(row)

# Material models from pfw_materials and schema docs
mat_text=read(pathlib.Path('scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_materials.py'))
models=[]
mg=re.search(r'MATERIAL_STRING_GENERATORS\s*=\s*\{(.*?)\}', mat_text, re.S)
if mg:
    models=re.findall(r"['\"]([^'\"]+)['\"]", mg.group(1))
# count presets by model
preset_models=[]
for m in re.finditer(r"([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*['\"]model['\"]\s*\]\s*=\s*['\"]([^'\"]+)['\"]", mat_text):
    preset_models.append((m.group(1),m.group(2)))
model_counts=Counter(model for _,model in preset_models)
# additional schema constitutive files relevant, all docs that have common MPM names or cohesive
schema_docs=repo/'src/coreComponents/schema/docs'
cohesive_files=[p.stem for p in schema_docs.glob('*Cohesive*.rst') if not p.stem.endswith('_other')]
# get schema attributes for selected material docs
material_docs=[]
selected=set(models)|set(cohesive_files)|{'ElasticIsotropic','ElasticTransverseIsotropicPressureDependent','CoupledCohesiveZone','BicrystalCohesiveZone','PolymerCohesiveZone'}
for name in sorted(selected):
    f=schema_docs/f'{name}.rst'
    if f.exists():
        text=f.read_text(errors='ignore')
        # parse rst grid table rows: name type default desc first columns (simple)
        attrs=[]
        for line in text.splitlines():
            if not line.strip() or set(line.strip())<=set('= '): continue
            parts=re.split(r'\s{2,}', line.strip(), maxsplit=3)
            if len(parts)>=3 and parts[0] not in ('Name',):
                attrs.append({'name':parts[0],'type':parts[1] if len(parts)>1 else '', 'default':parts[2] if len(parts)>2 else '', 'description':parts[3] if len(parts)>3 else ''})
        material_docs.append({'name':name,'preset_count':model_counts.get(name,0),'attributes':attrs[:12],'attribute_count':len(attrs),'source':f'src/coreComponents/schema/docs/{name}.rst'})
    else:
        material_docs.append({'name':name,'preset_count':model_counts.get(name,0),'attributes':[],'attribute_count':0,'source':''})

# Do not reuse cached MPM schema tables as current defaults. These references
# are derived from the same wrapper registrations used in the manual.
def write_wrapper_rst(name, rows):
    lines=[".. Generated from current source wrappers by extract_mpm_metadata.py.",
           ".. Full inherited XML schema generation requires a GEOS build.", "",
           "Current directly registered input attributes are listed below. Defaults",
           "are resolved from source where possible; remaining C++ expressions are",
           "shown verbatim. Inherited attributes are not repeated.", "",
           ".. list-table:: "+name+" source wrapper reference",
           "   :header-rows: 1", "   :widths: 24 16 18 42", "",
           "   * - Attribute", "     - Input / type", "     - Default", "     - Description"]
    for row in sorted(rows,key=lambda r:r["name"]):
        if row["flag"]=="FALSE":continue
        lines += ["   * - "+row["name"],
                  "     - "+(row["flag"]+" / "+row["type"]).strip(" /"),
                  "     - "+(row["default"] or ("required" if row["flag"]=="REQUIRED" else "not explicitly set")),
                  "     - "+row["description"].replace("\n"," ")]
    (schema_docs/(name+".rst")).write_text("\n".join(line.rstrip() for line in lines)+"\n")
if args.write_schema_docs:
    write_wrapper_rst("SolidMechanics_MPM", solver_wrappers)
    write_wrapper_rst("ParticleRegion", particle_region_wrappers)

# Examples/suite cases
suite_root=repo/'scripts/preProcessing/materialPointMethodParticleFileWriter'
def case_list(subdir):
    root=suite_root/subdir
    cases=[]
    for f in sorted(root.rglob('pfw_input_*.py')):
        rel=f.relative_to(root)
        name=str(rel.with_suffix('')).replace('pfw_input_','')
        family=str(rel.parent) if str(rel.parent)!='.' else ''
        # top comment block first 20 lines
        txt=f.read_text(errors='ignore')
        desc=''
        lines=[]
        for line in txt.splitlines()[:80]:
            if line.strip().startswith('#'):
                s=line.strip().lstrip('#').strip()
                if s: lines.append(s)
            elif lines:
                break
        desc=' '.join(lines[:3])[:240]
        cases.append({'family':family,'case':name,'path':str(f.relative_to(repo)),'description':desc})
    return cases
suite_cases={k:case_list(k) for k in ['examples','verification','validation']}

# Postprocessing scripts
post_scripts=[]
for f in sorted(suite_root.rglob('*')):
    if f.is_file() and (f.name.startswith('plot') or f.name.startswith('postProcess') or f.name.startswith('visitRender') or f.name in ['pfw_visit_render.py','mpm_vv_postprocess.py','mpm_vv_plot_tools.py','pfw_analysis.py','mpm_example_postprocess.py']):
        post_scripts.append(str(f.relative_to(repo)))
# More analysis scripts outside PFW? maybe benchmarks. Add scripts with postprocess/analyze/plot.
for f in sorted((repo/'scripts').glob('*.py')):
    if re.search(r'(plot|analysis|post|visit|paraview)', f.name, re.I):
        post_scripts.append(str(f.relative_to(repo)))
post_scripts=sorted(set(post_scripts))

# Source inventory
source_inventory=[
    {'component':'MPM solver','path':'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.[ch]pp'},
    {'component':'Dimension-specialized MPM mapping','path':'src/coreComponents/physicsSolvers/solidMechanics/MPMMappingKernels.hpp'},
    {'component':'Native interface measures','path':'src/coreComponents/physicsSolvers/solidMechanics/MPMGeometryMeasures.hpp'},
    {'component':'Axisymmetric ring measures','path':'src/coreComponents/mesh/AxisymmetricGeometry.hpp'},
    {'component':'MPM unified geometry','path':'src/coreComponents/mesh/ParticleGeometry.hpp'},
    {'component':'MPM particle subdivision','path':'src/coreComponents/mesh/ParticleSubdivision.hpp'},
    {'component':'MPM solver fields','path':'src/coreComponents/physicsSolvers/solidMechanics/MPMSolverFields.hpp'},
    {'component':'MPM events','path':'src/coreComponents/events/mpmEvents/*'},
    {'component':'Particle mesh and regions','path':'src/coreComponents/mesh/Particle*.{hpp,cpp}, src/coreComponents/mesh/particleGenerators/*'},
    {'component':'Solver enums','path':'src/coreComponents/physicsSolvers/solidMechanics/MPMSolverEnums.hpp'},
    {'component':'ParticleFileWriter','path':'scripts/preProcessing/materialPointMethodParticleFileWriter/particleFileWriter.py'},
    {'component':'PFW particle layout','path':'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_particle_layout.py'},
    {'component':'PFW geometry objects','path':'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_geometryObjects.py'},
    {'component':'PFW material library','path':'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_materials.py'},
    {'component':'PFW suite tools','path':'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_suite.py and mpm_vv_*.py'},
    {'component':'Canonical XML examples','path':'inputFiles/materialPointMethod/*'},
]

# Version
version=(repo/'src/VERSION').read_text(errors='ignore').strip() if (repo/'src/VERSION').exists() else ''

data={
 'version':version,
 'solver_wrappers':solver_wrappers,
 'event_data':event_data,
 'particle_types':particle_types,
 'particle_columns':particle_columns_data,
 'solver_enums':enums,
 'particle_region_wrappers':particle_region_wrappers,
 'pfw_params':params,
 'pfw_silo_common':silo_common,
 'pfw_silo_all':silo_all,
 'geometry_objects':geometry_names,
 'geometry_wrappers':wrapper_names,
 'geometry_utilities':utility_classes,
 'material_docs':material_docs,
 'material_preset_counts':dict(model_counts),
 'suite_cases':suite_cases,
 'post_scripts':post_scripts,
 'source_inventory':source_inventory,
}
fingerprint_paths=[
 'src/coreComponents/physicsSolvers/solidMechanics/MPMMappingKernels.hpp',
 'src/coreComponents/physicsSolvers/solidMechanics/MPMGeometryMeasures.hpp',
 'src/coreComponents/mesh/generators/InternalMeshGenerator.cpp',
 'src/coreComponents/mesh/mpiCommunications/SpatialPartition.cpp',
 'src/coreComponents/fileIO/silo/SiloFile.cpp',
 'src/coreComponents/fileIO/vtk/VTKPolyDataWriterInterface.cpp',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_geometryObjects.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_visit_render.py',
 'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.cpp',
 'src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp',
 'src/coreComponents/physicsSolvers/solidMechanics/MPMSolverEnums.hpp',
 'src/coreComponents/mesh/particleGenerators/ParticleMeshGenerator.cpp',
 'src/coreComponents/mesh/particleGenerators/ParticleMeshGenerator.hpp',
 'src/coreComponents/mesh/ParticleRegion.cpp',
 'src/coreComponents/mesh/ParticleGeometry.hpp',
 'src/coreComponents/mesh/AxisymmetricGeometry.hpp',
 'src/coreComponents/mesh/ParticleSubdivision.hpp',
 'src/coreComponents/mesh/ParticleSubRegionBase.hpp',
 'src/coreComponents/mesh/ParticleSubRegionBase.cpp',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/particleFileWriter.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_surface_refinement.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_particle_layout.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_simplex_clipping.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_reference_geometry.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_axisymmetric.py',
 'scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_check.py',
]
data['source_sha256']={p:hashlib.sha256((repo/p).read_bytes()).hexdigest() for p in fingerprint_paths}
path=args.out
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(data,indent=2), encoding='utf-8')
print('wrote', path)
print('solver wrappers',len(solver_wrappers),'events',len(event_data),'pfw params',len(params),'geometry',len(geometry_names),'wrappers',len(wrapper_names),'materials',len(material_docs),'post scripts',len(post_scripts))
print('event names:', ', '.join(e['name'] for e in event_data))
