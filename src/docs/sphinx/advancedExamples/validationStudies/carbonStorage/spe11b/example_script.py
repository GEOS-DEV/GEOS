#SPHINX_IMPORT_BEGIN
import os, logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#SPHINX_IMPORT_END

#SPHINX_DATA_BEGIN
# for geos
name_indirection = {'pressure': 'pres',
                         'phaseVolumeFraction_0': 'satg',
                         'phaseVolumeFraction_1': 'satw',
                         'fluid_phaseCompFraction_2': 'mCO2',
                         'fluid_phaseCompFraction_1': 'mH2O',
                         'fluid_phaseDensity_0': 'rG',
                         'fluid_phaseDensity_1': 'rL',
                         'rockPorosity_porosity': 'poro',
                         'elementVolume': 'vol',
                         'phaseMobility_0': 'krg',
                         'phaseMobility_1': 'krw'}

logging.basicConfig(filename='spe11-pt.log', encoding='utf-8', level=logging.INFO)


SEC2YEAR =  365*24*3600

#main containers
class Data:

    data_sets = {}
    schedule = []
    PO1 = []
    PO2 = []
    boxes = {}

data = Data()
def define_boxes():
    #define and rescale
    data.boxes = {'Whole': [(0.0, 0.0, 0.0), (8400, 1., 1200.)]}

    data.boxes['A'] = [(3300.,0.,0.), (8300.,1.,600.)]
    data.boxes['B'] = [(100.,0.,600.), (3300.,1.,1200.)]
    data.PO1 = [ 4500, 500]
    data.PO2 = [ 5100, 1100]
    data.schedule = np.arange(0., 1000 * SEC2YEAR, 1000/200 * SEC2YEAR)
#SPHINX_DATA_END

#SPHINX_UTILS_BEGIN
formula = {'mImmobile': 'if(krg<0.001, rG*satg*poro*vol)',
           'mMobile': 'if(krg>0.001,rG*satg*poro*vol)',
           'mDissolved': 'rL*mCO2*poro*vol*satw' }
def process_keys(formula, fielddict):
    """
    Basic formula interface, reading from multiplication and division (no parenthesis)
    :param formula: string of the formula
    :param fielddict: all-containing dict with base-variables
    :return: Array of the formula value
    """
    import re
    if re.findall(r'if', formula):
        bits = formula.split(',')
        assert (len(bits) == 2)
        cdt = bits[0][3:]
        if re.findall(r'>', cdt):
            inner_bits = cdt.split('>')
            cdt_field = inner_bits[0].strip()
            cdt_lim = float(inner_bits[1])
            return np.where(fielddict[cdt_field] > cdt_lim, process_keys(bits[1][:-1], fielddict), 0)
        elif re.findall(r'<', cdt):
            inner_bits = cdt.split('<')
            cdt_field = inner_bits[0].strip()
            cdt_lim = float(inner_bits[1])
            return np.where(fielddict[cdt_field] < cdt_lim, process_keys(bits[1][:-1], fielddict), 0)
        else:
            raise NotImplemented('Not a valid condition')
    elif re.findall(r'\+', formula):
        bits = formula.split('+')
        bits = [item.strip() for item in bits]
        return process_keys(bits[0], fielddict) + process_keys('+'.join(bits[1:]), fielddict)
    elif re.findall(r'-', formula):
        bits = formula.split('-')
        bits = [item.strip() for item in bits]
        return process_keys(bits[0], fielddict) - process_keys('+'.join(bits[1:]), fielddict)
    elif re.findall(r'\*', formula):
        bits = formula.split('*')
        bits = [item.strip() for item in bits]
        return process_keys(bits[0], fielddict) * process_keys('*'.join(bits[1:]), fielddict)
    elif re.findall(r'/', formula):
        bits = formula.split('/')
        bits = [item.strip() for item in bits]
        return process_keys(bits[0], fielddict) / process_keys('/'.join(bits[1:]), fielddict)
    return fielddict[formula]

#SPHINX_UTILS_END

#SPHINX_INTERPOLATE_BEGIN
def get_interpolate(points_from_vtk, fields: dict, nskip=1):
    """ getting dict of proper interpolation for fields """
    def get_lambda_2(xpts, disc_func):
        from scipy.interpolate import LinearNDInterpolator
        x, z = xpts
        return LinearNDInterpolator(np.asarray([x, z]).transpose(), disc_func, fill_value=0.0)

    fn = {}
    for key, val in fields.items():
        if key in name_indirection.values() or key in formula.keys():
            logging.info(f'Building interpolation for {key}')
            # do not interpolate if not interested by (and only brought as field is another component of
            # a field of interest
            T = points_from_vtk
            fn[key] = get_lambda_2((T[::nskip, 0], T[::nskip, 2]), val[::nskip])

    return fn

def integrate_2(pts, fields, box):
    """ Simplistic integration equivalent to Paraview integrate variables """
    xmin, ymin, zmin = box[0]
    xmax, ymax, zmax = box[1]
    ii = np.intersect1d(np.intersect1d(np.argwhere((pts[:, 0] > xmin)), np.argwhere((pts[:, 0] < xmax))),
                        np.intersect1d(np.argwhere((pts[:, 2] > zmin)), np.argwhere((pts[:, 2] < zmax))))
    return fields[ii].sum()

def intergrate_fields(pft_tuple: list):

    logging.info(f'Interpolating reported fields')
    lines = list()
    pts_from_vtk, fields, time = pft_tuple

    # some lines for MC magic number
    for key, form in formula.items():
        fields[key] = process_keys(form, fields)

    fn = get_interpolate(pts_from_vtk,
                                    {'pres': fields['pres'], 'vol': fields['vol']},
                                    nskip=1)

    line = [time]
    logging.info(f'Writing sparse time {time:2}')
    # deal with P1 and P2
    p1 = fn['pres'](data.PO1)
    p2 = fn['pres'](data.PO2)
    line.extend([p1[0], p2[0]])
    # deal box A & B
    logging.info(f'Integrating for boxes A and B')
    for box_name, box in data.boxes.items():
            if box_name in ['A', 'B']:
                line.extend([
                    integrate_2(pts_from_vtk, fields['mMobile'], box),
                    integrate_2(pts_from_vtk, fields['mImmobile'], box),
                    integrate_2(pts_from_vtk, fields['mDissolved'], box)
                ])
        # #deal sealTot
    lines.append(line)

    return pd.DataFrame(data=lines,
                        columns=['t[s]', 'p1[Pa]', 'p2[Pa]', 'mobA[kg]', 'immA[kg]', 'dissA[kg]',
                                 'mobB[kg]', 'immB[kg]', 'dissB[kg]'])
#SPHINX_INTERPOLATE_END

#SPHINX_PTIME_BEGIN
def process_time(pvdfile, time, olist):

    import vtk
    import numpy as np
    from vtk.util import numpy_support

    reader = vtk.vtkXMLMultiBlockDataReader()
    reader.SetFileName( './' + os.path.dirname(pvdfile) + '/' + data.data_sets[time] )
    loadlist = list(olist)

    for fieldname in loadlist:
        if fieldname in olist or fieldname in ['ghostRank']:
            reader.SetCellArrayStatus(fieldname, 1)
        else:
            reader.SetCellArrayStatus(fieldname, 0)

    logging.info('Updating multiblocks')
    reader.Update()
    logging.info('Processing multiblocks')
    # remove well blocks - in hierarchy mesh.lvl0.CellElementRegion/WellElementRegion
    if reader.GetOutput().GetBlock(0).GetBlock(0).GetNumberOfBlocks() > 1:
        reader.GetOutput().GetBlock(0).GetBlock(0).RemoveBlock(1)

    # browse blocks - assemble size
    # counting cells and blocks/rank
    nv = reader.GetOutput().GetNumberOfCells()
    it = reader.GetOutput().NewIterator()
    it.InitTraversal()

    ncnf = 0
    it = reader.GetOutput().NewIterator()
    it.InitTraversal()
    for ifields in olist:
        ncc = reader.GetOutput().GetDataSet(it).GetCellData().GetArray(ifields).GetNumberOfComponents()
        ncnf += ncc
        it.GoToNextItem()

    # form data container
    f = np.zeros(shape=(nv, ncnf), dtype='float')
    pts = np.zeros(shape=(nv, 3), dtype='float')

    # get seal flagged
    mesh = reader.GetOutput().GetBlock(0).GetBlock(0).GetBlock(0)
    fielddict = {}
    for i in range(0, mesh.GetNumberOfBlocks()):
        start = 0
        it = mesh.GetBlock(i).NewIterator()
        it.InitTraversal()
        info = mesh.GetMetaData(i)
        while not it.IsDoneWithTraversal():
            field = mesh.GetBlock(i).GetDataSet(it).GetCellData().GetArray("ghostRank")
            it.GoToNextItem()

        nc = 0
        j = 0  # for field loop
        for ifields in olist:
            start = 0
            it = reader.GetOutput().NewIterator()
            it.InitTraversal()
            while not it.IsDoneWithTraversal():
                field = reader.GetOutput().GetDataSet(it).GetCellData().GetArray(ifields)
                nt = field.GetNumberOfValues()
                nc = field.GetNumberOfComponents()
                f[start:(start + int(nt / nc)), j:j + nc] = numpy_support.vtk_to_numpy(field).reshape(int(nt / nc),
                                                                                                      nc)

                cc = vtk.vtkCellCenters()
                cc.SetInputData(reader.GetOutput().GetDataSet(it))
                cc.Update()
                if nt > 0:
                    pts[start:(start + int(nt / nc)), :] = numpy_support.vtk_to_numpy(
                        cc.GetOutput().GetPoints().GetData())
                it.GoToNextItem()
                start += int(nt / nc)

            if nc > 1:
                for k in range(0, nc):
                    if ifields + "_" + str(k) in name_indirection:
                        fielddict[name_indirection[ifields + "_" + str(k)]] = f[:, j + k]
            else:
                if ifields in name_indirection:
                    fielddict[name_indirection[ifields]] = f[:, j]
            j = j + nc

    return pts, fielddict, time
#SPHINX_PTIME_END



#SPHINX_MAIN_BEGIN
def write(ifile):

    #utils
    def read_pvd(ifile):

        logging.info(f'Reading pvd file {ifile}')
        import xml.etree.ElementTree as ET
        tree = ET.parse(ifile)
        root = tree.getroot()
        for ds in root.find('Collection').findall('DataSet'):
            data.data_sets[float(ds.attrib['timestep'])] = ds.attrib['file']


    # note that solubility is in KgCO2/KgBrine
    read_pvd(ifile)
    define_boxes()
    # preprocess input list from desired output
    import re
    olist_ = set(
        [item for name in name_indirection.keys() for item in re.findall(r'(\w+)_\d|(\w+)', name)[0] if
         len(item) > 0])

    ddf = pd.DataFrame()
    for time in data.schedule:
        field_and_pts = process_time(ifile, time, olist_)
        ddf = pd.concat([ddf, intergrate_fields(field_and_pts)])
    ddf.to_csv( './spe11b_time_series.csv',
                       index_label=False)



def plot():
    df = pd.read_csv('./spe11b_time_series.csv')
    fig, axs = plt.subplots(2, 2)
    (time_name, time_unit), (mass_name, mass_unit), (pressure_name, pressure_unit) \
        = ('Y', SEC2YEAR), ('T', 1000), ('bar', 1e5)
    # pressures
    axs[0][0].plot(df['t[s]'].to_numpy() / time_unit, df['p1[Pa]'].to_numpy() / pressure_unit,
                   label=f'pressure 1 [{pressure_name}]')
    axs[0][0].plot(df['t[s]'].to_numpy() / time_unit, df['p2[Pa]'].to_numpy() / pressure_unit,
                   label=f'pressure 2 [{pressure_name}]')
    axs[0][0].legend()
    # box A
    axs[0][1].plot(df['t[s]'].to_numpy() / time_unit, df['mobA[kg]'].to_numpy() / mass_unit,
                   label=f'mobile CO2 [{mass_name}]')
    axs[0][1].plot(df['t[s]'].to_numpy() / time_unit, df['immA[kg]'].to_numpy() / mass_unit,
                   label=f'immobile CO2 [{mass_name}]')
    axs[0][1].plot(df['t[s]'].to_numpy() / time_unit, df['dissA[kg]'].to_numpy() / mass_unit,
                   label=f'dissolved CO2 [{mass_name}]')
    axs[0][1].legend()
    axs[0][1].set_title('boxA')
    # box B
    axs[1][0].plot(df['t[s]'].to_numpy() / time_unit, df['mobB[kg]'].to_numpy() / mass_unit,
                   label=f'mobile CO2 [{mass_name}]')
    axs[1][0].plot(df['t[s]'].to_numpy() / time_unit, df['immB[kg]'].to_numpy() / mass_unit,
                   label=f'immobile CO2 [{mass_name}]')
    axs[1][0].plot(df['t[s]'].to_numpy() / time_unit, df['dissB[kg]'].to_numpy() / mass_unit,
                   label=f'dissolved CO2 [{mass_name}]')
    axs[1][0].legend()
    axs[1][0].set_title('boxB')

    fig.tight_layout()
    fig.savefig('./spe11b_timeseries.png', bbox_inches='tight')


if __name__ == "__main__":

    write('vtkOutput.pvd')
    # plot
    plot()
#SPHINX_MAIN_END
