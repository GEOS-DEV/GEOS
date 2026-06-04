import re, json, pathlib, ast
from collections import defaultdict, Counter
repo=pathlib.Path('/mnt/data/geos_mpm_review/GEOS-mpm')

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
        rows.append({'name':name,'flag':flag,'default':default,'restart':restart,'description':desc,'source':str(path),'line':line_no(text,idx)})
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
# defaults for columns from source, manually parse switch cases somewhat
# Use hard-coded defaults based on observed switch cases
col_defaults={c:'required' for c in particle_columns}
for c in ['StrengthScale','MaterialDirectionXX','MaterialDirectionYY','MaterialDirectionZZ']:
    col_defaults[c]='1.0 if omitted'
for c in ['Temperature']:
    col_defaults[c]='300.0 if omitted'
for c in ['ParticleType']:
    col_defaults[c]='2.0 (CPDI) if omitted'
for c in ['CZTag','MaterialType','ContactGroup','Damage','Porosity','VelocityX','VelocityY','VelocityZ','MaterialDirectionXY','MaterialDirectionXZ','MaterialDirectionYX','MaterialDirectionYZ','MaterialDirectionZX','MaterialDirectionZY','SurfaceNormalX','SurfaceNormalY','SurfaceNormalZ','SurfacePositionX','SurfacePositionY','SurfacePositionZ','SurfaceTractionX','SurfaceTractionY','SurfaceTractionZ','TemperatureRate']:
    if c in col_defaults: col_defaults[c]='0.0 if omitted'
particle_columns_data=[{'name':c,'default':col_defaults[c]} for c in particle_columns]

# Solver enum options from SolidMechanicsMPM.hpp
htext=read(pathlib.Path('src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp'))
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
    enums.append({'name':name,'values':vals})

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

# Schema docs for ParticleMesh, ParticleRegion, SolidMechanics_MPM parse
schema_tables={}
for name in ['SolidMechanics_MPM','ParticleMesh','ParticleRegion','MPMEvents']:
    f=schema_docs/f'{name}.rst'
    if f.exists():
        rows=[]
        for line in f.read_text(errors='ignore').splitlines():
            if not line.strip() or set(line.strip())<=set('= '): continue
            parts=re.split(r'\s{2,}', line.strip(), maxsplit=3)
            if len(parts)>=3 and parts[0]!='Name':
                rows.append({'name':parts[0], 'type':parts[1], 'default':parts[2], 'description':parts[3] if len(parts)>3 else ''})
        schema_tables[name]=rows

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
    {'component':'MPM kernels','path':'src/coreComponents/physicsSolvers/solidMechanics/kernels/ExplicitMPM.hpp'},
    {'component':'MPM solver fields','path':'src/coreComponents/physicsSolvers/solidMechanics/MPMSolverFields.hpp'},
    {'component':'MPM events','path':'src/coreComponents/events/mpmEvents/*'},
    {'component':'Particle mesh and regions','path':'src/coreComponents/mesh/Particle*.{hpp,cpp}, src/coreComponents/mesh/particleGenerators/*'},
    {'component':'Schema-generated docs','path':'src/coreComponents/schema/docs/SolidMechanics_MPM.rst and related MPM docs'},
    {'component':'ParticleFileWriter','path':'scripts/preProcessing/materialPointMethodParticleFileWriter/particleFileWriter.py'},
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
 'pfw_params':params,
 'pfw_silo_common':silo_common,
 'pfw_silo_all':silo_all,
 'geometry_objects':geometry_names,
 'geometry_wrappers':wrapper_names,
 'geometry_utilities':utility_classes,
 'material_docs':material_docs,
 'material_preset_counts':dict(model_counts),
 'schema_tables':schema_tables,
 'suite_cases':suite_cases,
 'post_scripts':post_scripts,
 'source_inventory':source_inventory,
}
path=pathlib.Path('/mnt/data/geos_mpm_extracted.json')
path.write_text(json.dumps(data,indent=2))
print('wrote', path)
print('solver wrappers',len(solver_wrappers),'events',len(event_data),'pfw params',len(params),'geometry',len(geometry_names),'wrappers',len(wrapper_names),'materials',len(material_docs),'post scripts',len(post_scripts))
print('event names:', ', '.join(e['name'] for e in event_data))
