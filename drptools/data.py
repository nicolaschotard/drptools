from __future__ import print_function
import os
import glob
import numpy
import fitsio
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, Column
from astropy.units import Quantity
from termcolor import colored
from . import utils as cutils

try:
    from lsst.afw import image as afwimage
    from lsst.afw import table as afwtable
    import lsst.daf.persistence as dafPersist
except ImportError:
    print(colored("WARNING: LSST stack is probably not installed", "yellow"))


class DRPLoader(object):

    """Load an LSST DRP output and a few useful things."""
    
    def __init__(self, drp_path):
        """The 'drp_path' input is the path to the DRP output directory."""

        # Load the bulter
        self.drp_path = drp_path
        self.butler = dafPersist.Butler(self.drp_path)
        self.mapper = self.butler._getDefaultMapper()
        self.repoData = self.butler._repos.outputs()[0]

        # Load some basic info on the current DRP
        self.repo_input = self._get_repo("input")    # repo containing the raw data after ingestion
        self.repo_output = self._get_repo("output")  # repo containing the processed data

        # Load some dataids
        self.datasetTypes = self._get_datasetTypes()
        self.datasetTypes_filename = self._get_datasetTypes_withfiles()
        self.dataIds = {}
        for dataset in ['raw', 'forced_src', 'deepCoadd_meas', 'deepCoadd_forced_src']:
            dataids = self.get_dataIds(dataset)
            if len(dataids):
                self.dataIds[dataset] = dataids

        # Load filter and visits
        self.filters = self.get_filter_list()
        self.visits = self.get_visit_list()

        # Skymap if any
        if self.butler._locate('deepCoadd_skyMap', dafPersist.DataId({}), write=False) is not None:
            self.skymap = self.butler.get("deepCoadd_skyMap")
            self.skymap_name = self.skymap.__class__.__name__
            self.skymap_doc = self.skymap.__doc__
            self.skymap_config = self.skymap.config.toDict()
            self.skymap_numtracts = self.skymap._numTracts
            self.skymap_numpatches = self.skymap[0].getNumPatches()

        # Mapper info
        self.mapper_name = self.mapper.__name__
        self.mapper_package = self.mapper.packageName 
        self.mapper_camera = self.mapper.getCameraName()

        # Packages
        self.packages = self.butler.get('packages')

        # Other
        self.configs = self._load_configs()
        self.schemas = self._load_schemas()

    def _get_repo(self, repo):
        """Get the full path of the input/output repository."""
        has_parent = bool(len(self.repoData.getParentRepoDatas()))
        if repo == 'output':
            # The given repo is an output one if it has a parent, otherwise, it should be an input one
            return self.repoData.cfgRoot if has_parent else 'No output directory found'
        elif repo == 'input':
            # the parent directory is the one containing the raw data (after ingestion)
            if has_parent:
                parentRepoData = self.repoData.getParentRepoDatas()[0]         
                return os.path.realpath(parentRepoData.cfgRoot)  # input path -> ../input in this case
            else:
                return self.repoData.cfgRoot
        else:
            raise IOError("Wrong repo name. You should not be calling this internal method anyway.")

    def _get_datasetTypes(self):
        return sorted(self.repoData.repo._mapper.mappings.keys())

    def _get_datasetTypes_withfiles(self):
        mappables = [m for m in dir(self.repoData.repo._mapper) if m.startswith('map_')]
        withfile = [m.replace('map_', '') for m in mappables if m.endswith('_filename')]
        return sorted(withfile)

    def get_catIdKeys(self, datasetType):
        """Get the list of ID keys for a given catalog."""
        if datasetType not in self.datasetTypes:
            msg = "%s is not a valid datasetType. Check self.datasetTypes for the valid list." % \
                  datasetType
            raise IOError(msg)
        return self.butler.getKeys(datasetType)

    def get_dataIds(self, datasetType):
        """Get all available data id for a given dataType."""
        keys = self.get_catIdKeys(datasetType)
        # 'tract' is present in the keys for the forced_src cat, but should not be
        if datasetType == 'forced_src':
            del keys['tract']
        try:
            metadata = self.butler.queryMetadata(datasetType, format=sorted(keys.keys()))
        except:
            metadata = None
        if metadata is not None:
            return [dict(zip(sorted(keys.keys()), list(v) if not isinstance(v, list) else v))
                    for v in metadata]
        else:
            if datasetType not in self.repoData.repo._mapper.datasets:
                return []
            template = self.repoData.repo._mapper.datasets[datasetType]._template
            path = os.path.join(self.repoData.cfgRoot, os.path.dirname(template))
            basepath = "/".join([p for p in path.split('/') if not p.startswith('%')]) + "/"
            keys = self.butler.getKeys(datasetType)
            gpath = "/".join([p if not p.startswith('%') else '*' for p in path.split('/')])
            paths = [p for p in glob.glob(gpath) if 'merged' not in p]
            return [{k: keys[k](v) for k, v in zip(keys, p.split(basepath)[1].split('/'))}
                    for p in paths]

    def get_filter_list(self):
        """Get the list of filters."""
        return set([dataid['filter'] for dataid in self.dataIds['raw']])

    def get_visit_list(self):
        """"All available vists."""
        visits = {filt: list(set([dataid['visit'] 
                                  for dataid in self.dataIds['raw'] if dataid['filter'] == filt]))
                  for filt in self.filters}
        return visits

    def _load_configs(self):
        """Load configs for the main tasks."""
        configs = self._load_generic_dataset("config")
        return {cfg: configs[cfg].toDict() for cfg in configs}

    def _load_schemas(self):
        """Load the schemas for all catalogs."""
        schemas = self._load_generic_dataset("schema")
        for schema in schemas:
            sch = schemas[schema].asAstropy()
            schemas[schema] = {col: {'description': sch[col].description, 
                                     'unit': sch[col].unit,
                                     'dtype': sch[col].dtype}
                               for col in sch.colnames}
        return schemas
        
    def _load_generic_dataset(self, datatype):
        """Load the schema or config datasets."""
        if datatype not in ['config', 'schema']:
            raise IOError("`datatype` must be either `config` or `schema`.")
        datatypes = {}
        for dataset in self.datasetTypes:
            if not dataset.endswith('_%s' % datatype):
                continue
            for dataId in ([{}] + [self.dataIds[key][0] for key in self.dataIds]):
                try:
                    datatypes[dataset] = self.butler.get(dataset, dataId=dataId)
                    break
                except:
                    pass
        return datatypes
    
    def overview(self):
        """Overview of the current DRP content."""
        # General info on the data repository
        html = "<h2>General info on the current DRP</h2>"
        
        # Repos
        html += "<h3>Paths to the repositories</h3>")
        html += cutils.square_list(["<b>Input</b>: %s</li>" % BUTLER.repo_input,
                                    "<b>Output</b>: %s</li>" % BUTLER.repo_output
                                   ]
                                  )

        # Info on the mapper, camera, package
        html += "<h3>Mapper info</h3>"
        html += cutils.square_list(["<b>Package</b>: %s" % BUTLER.mapper_package,
                                     "<b>Camera</b>: %s" % BUTLER.mapper_camera,
                                     "<b>Name</b>: %s" % BUTLER.mapper_name
                                    ]
                                   )
                                   
            
        html += "<h3>Filters and visits</h3>"
        html += "<table>"
        html += "<tr><th>Name</th>"
        html += "".join(["<td>%s</td>" % filt for filt in BUTLER.filters])
        html += "</tr>"
        html += "<tr><th>#Visits</th>"
        html += "".join(["<td>%i</td>" % len(BUTLER.visits[filt]) for filt in BUTLER.filters])
        html += "</tr>"
        html += "</table>"
        
        # Other info, filter, skymap, etc.
        items = []
        if hasattr(self, 'skymap'):
            items.append("<b>Sky map</b>: %s" % str(self.skymap_name))
        if len(items):
            html += "<h3>Other info</h3>"
            html += cutils.square_list(items)

        return page


class DRPFiles(DRPLoader):

    def __init__(self, drp_path):
        """The 'drp_path' input is the path to the DRP output directory."""
        super().__init__(drp_path)

    def get_file(self, datatype, dataid):
        try:
            cfiles = self.butler.get('%s_filename' % datatype, dataId=dataid)
            for i, cfile in enumerate(cfiles):
                if self.repo_output in cfile and not os.path.exists(cfile):
                    cfiles[i] = cfile.replace(self.repo_output, self.repo_input)
            return cfiles
        except:
            return []

    def get_files(self, datatype, filt=None, visit=None, tract=None, patch=None):
        dataids = self.get_dataid_from_dataset(datatype)
        files = numpy.concatenate([self.get_file(datatype, dataid) for dataid in dataids])
        return files

    def get_dataid_from_dataset(self, datatype, test=False):
        try:
            keys = self.butler.getKeys(datatype)
        except:
            keys = {}
        if not len(keys):
            return [{}]
        elif 'visit' in keys and 'tract' in keys:
            key = 'forced_src'
            dataIds = [self.dataIds[key][0]] if test else self.dataIds[key]
        elif 'visit' in keys:
            key = 'raw'
            dataIds = [self.dataIds[key][0]] if test else self.dataIds[key]
        elif 'tract' in keys:
            key = 'deepCoadd_meas'
            dataIds = [self.dataIds[key][0]] if test else self.dataIds[key]
        else:
            dataIds = self.get_dataIds(datatype)
            if not len(dataIds):
                dataIds = [{}]
        return [{k: dataid[k] for k in dataid if k in keys} for dataid in dataIds]

    def _has_file(self, datatype):
        return bool(len(self.get_file(datatype, self.get_dataid_from_dataset(datatype, test=True)[0])))

    def _type_file(self, datatype):
        return self.get_file(datatype, self.get_dataid_from_dataset(datatype, test=True)[0])


class DRPCatalogs(DRPLoader):

    """Load catalogs from an LSST DRP output path."""

    def __init__(self, drp_path):
        """The 'drp_path' input is the path to the DRP output directory."""
        super().__init__(drp_path)

        # Initialize data dictionnaries
        self.dataids = {}
        self.catalogs = {}
        self.keys = {}
        self.missing = {}
        self.from_butler = {'getmag': None, 'wcs': None,
                            'schema': None, 'extension': None}
        self.append = False

    def _load_dataids(self, catalog, **kwargs):
        """Get the 'forced_src' catalogs."""
        print("INFO: Getting list of available data for", catalog)
        if 'deepCoadd' in catalog:  # The deepCoadd* catalogs
            deepcoadd = [cat for cat in self.dataids if 'deepCoadd' in cat]
            if len(deepcoadd):
                dataids = self.dataids[deepcoadd[0]]
            else:
                dataids = [dict(tract=tract.getId(), patch="%d,%d" % patch.getIndex(), filter=filt)
                           for tract in self.butler.get("deepCoadd_skyMap")
                           for patch in tract
                           for filt in kwargs.get('filter', ['u', 'g', 'r', 'i', 'i2', 'z'])]
        else:  # The other catalogs
            keys = self.butler.getKeys(catalog)
            if 'tract' in keys:
                keys.pop('tract')
                metadata = self.butler.queryMetadata(
                    catalog, format=sorted(keys.keys()))
                dataids = [cutils.merge_dicts(dict(zip(sorted(keys.keys()), list(v))), {'tract': 0})
                           for v in metadata]
            else:
                metadata = self.butler.queryMetadata(
                    catalog, format=sorted(keys.keys()))
                dataids = [dict(zip(sorted(keys.keys()), [v] if not isinstance(v, list) else v))
                           for v in metadata]

        if len(dataids) == 0:
            raise IOError(
                "No dataIds. Check the catalog, the config file, and path to the bulter.")

        # Specific selection make by the user?
        for kwarg in kwargs:
            if kwarg not in dataids[0]:
                continue
            print("INFO: Selecting data ids according to the '%s' selection" % kwarg)
            print("  - input: %i data ids" % len(dataids))
            if not isinstance(kwargs[kwarg], list):
                kwargs[kwarg] = [kwargs[kwarg]]
            dataids = [dataid for dataid in dataids if dataid[kwarg]
                       in kwargs[kwarg]]
            print("  - selected: %i data ids" % len(dataids))

        # Select the ccd/visit according to the input list of patch if given
        if 'deepCoadd' not in catalog and 'patch' in kwargs and 'filter' in kwargs:
            print("INFO: Selecting visit/ccd according to the input list of patches")
            print("  - input: %i data ids" % len(dataids))
            ccds_visits = self._get_ccd_visits(**kwargs)
            dataids = [dataid for dataid in dataids if
                       (dataid['ccd'], dataid['visit']) in ccds_visits]
            print("  - selected: %i data ids" % len(dataids))

        # Only keep dataids with data
        print("INFO: Keep data IDs with data on disk")
        print("  - input: %i data ids" % len(dataids))
        self.dataids[catalog] = [dataid for dataid in dataids if
                                 self.butler.datasetExists(catalog, dataId=dataid)]
        self.missing[catalog] = [dataid for dataid in dataids if not
                                 self.butler.datasetExists(catalog, dataId=dataid)]
        print("  - selected: %i data ids" % len(self.dataids[catalog]))
        if len(self.missing[catalog]):
            print("  - missing: %i data ids (list available in 'self.missing[catalog]':" %
                  len(self.missing[catalog]))
        print("INFO: %i data ids finally kept" % len(self.dataids[catalog]))
        if len(self.dataids[catalog]) == 0:
            raise IOError(
                "No data found for this catalog. Remove this catalog from the list.")

    def _get_ccd_visits(self, **kwargs):
        """Return the available ccd/visit according to the input list of patch."""
        dids = [{'filter': filt, 'patch': patch, 'tract': 0}
                for filt in kwargs['filter'] for patch in kwargs['patch']
                if self.butler.datasetExists('deepCoadd',
                                             dataId={'filter': filt, 'patch': patch, 'tract': 0})]
        filenames = [kwargs['butler'] + '/deepCoadd' + "/%s/%i/%s.fits" %
                     (did['filter'], did['tract'], did['patch']) for did in dids]
        if self.from_butler['extension'] is None:
            # Extensions have no names, so we have to guess which extension contains the
            # 'visit' and 'ccd' info. we do it only once and then store this info. This allows
            # us to be safe against variations in the number of extensions in fits files.
            fitsdata = fitsio.FITS(filenames[0])
            self.from_butler['extension'] = [i for i, ext in enumerate(fitsdata)
                                             if (all([(key in ext.get_colnames())
                                                      for key in ['ccd', 'visit']])
                                                 if ext.get_exttype() == 'BINARY_TBL'
                                                 else False)][0]
            # other idea from Jim
            # fitsdata = fitsio.read(filenames[0], 4)
            # fitsdata[numpy.array([n.startswith('CoaddInputs')
            #                    for n in fitsdata['name']])]['cat.archive'][1] + 4
        return numpy.concatenate([fitsio.read(filename, columns=['ccd', 'visit'],
                                           ext=self.from_butler['extension'])
                               for filename in filenames]).tolist()

    def _load_catalog_dataid(self, catalog, dataid, table=True, **kwargs):
        """Load a catalog from a 'dataId' set of parameter."""
        try:
            cat = self.butler.get(catalog, dataId=dataid,
                                  flags=afwtable.SOURCE_IO_NO_FOOTPRINTS)
        except:  # OperationalError: no such column: flags
            cat = self.butler.get(catalog, dataId=dataid)
        if self.from_butler['schema'] is None and hasattr(cat, 'getSchema'):
            self.from_butler['schema'] = cat.getSchema()
        return cat.getColumnView().extract(*self.keys[catalog],
                                           copy=True, ordered=True) if table else cat

    def _get_catalog(self, dataset, **kwargs):
        """Load the catalogs from the butler."""
        filenames = (self.butler.get(dataset + "_filename",
                                     dataId, immediate=True)[0]
                     for dataId in self.dataids[dataset])
        try:  # In recent stack version, metadata are in HDU 1
            headers = (afwimage.readMetadata(fn, 1) for fn in filenames)
            size = sum(md.get("NAXIS2") for md in headers)
        except:  # Older stack version
            headers = (afwimage.readMetadata(fn, 2) for fn in filenames)
            size = sum(md.get("NAXIS2") for md in headers)
        cat = self.butler.get(dataset, self.dataids[dataset][0],
                              flags=afwtable.SOURCE_IO_NO_FOOTPRINTS, immediate=True)
        self.from_butler['schema'] = cat.schema
        catadic = {k: [] for k in sorted(self.dataids[dataset][0].keys())}
        catalog = afwtable.SourceCatalog(self.from_butler['schema'])
        catalog.reserve(size)
        pbar = cutils.progressbar(len(self.dataids[dataset]))
        print("INFO: Looping over the dataids")
        for i, dataid in enumerate(self.dataids[dataset]):
            cat = self.butler.get(dataset, dataid,
                                  flags=afwtable.SOURCE_IO_NO_FOOTPRINTS)
            catalog.extend(cat, deep=True)
            for newkey in catadic:
                catadic[newkey].extend([dataid[newkey]] * len(cat))
            pbar.update(i + 1)
        pbar.finish()
        print("INFO: Merging the dictionnaries")
        catadic.update(catalog.getColumnView().extract(*self.keys[dataset],
                                                       copy=True, ordered=True))
        # Clean memory before going further
        # gc.collect()
        return catadic

    def _load_catalog(self, catalog, **kwargs):
        """Load a given catalog."""
        self._load_dataids(catalog, **kwargs)
        print("INFO: Getting the data from the butler for %i fits files" %
              len(self.dataids[catalog]))
        self.catalogs[catalog] = Table(self._get_catalog(catalog, **kwargs))
        print("INFO: Getting descriptions and units")
        for k in self.catalogs[catalog].keys():
            if k in self.from_butler['schema']:
                asfield = self.from_butler['schema'][k].asField()
                self.catalogs[catalog][k].description = cutils.shorten(asfield.getDoc())
                self.catalogs[catalog][k].unit = asfield.getUnits()
        self.from_butler['schema'] = None
        print("INFO: %s catalog loaded (%i sources)" %
              (catalog, len(self.catalogs[catalog])))
        self._add_new_columns(catalog)
        if 'matchid' in kwargs and catalog == 'forced_src':
            self._match_ids()
        if 'output_name' in kwargs:
            self.save_catalogs(kwargs['output_name'], catalog,
                               kwargs.get('overwrite', False), delete_catalog=True)

    def _match_deepcoadd_catalogs(self):
        """In case of missing data for one catalog, remove corresonding data from the other."""
        if 'deepCoadd_meas' in self.catalogs and 'deepCoadd_forced_src' in self.catalogs:
            if len(self.catalogs['deepCoadd_meas']) == len(self.catalogs['deepCoadd_forced_src']):
                return
            print(colored("\nINFO: matching 'deepCoadd_meas' and 'deepCoadd_forced_src' catalogs",
                          'green'))
            for dataid in self.missing['deepCoadd_meas']:
                filt = (self.catalogs['deepCoadd_forced_src']['filter'] == dataid['filter']) & \
                       (self.catalogs['deepCoadd_forced_src']
                        ['patch'] == dataid['patch'])
                self.catalogs['deepCoadd_forced_src'] = self.catalogs['deepCoadd_forced_src'][~filt]
            for dataid in self.missing['deepCoadd_forced_src']:
                filt = (self.catalogs['deepCoadd_meas']['filter'] == dataid['filter']) & \
                       (self.catalogs['deepCoadd_meas']
                        ['patch'] == dataid['patch'])
                self.catalogs['deepCoadd_meas'] = self.catalogs['deepCoadd_meas'][~filt]

    def _match_ids(self):
        """Select in the 'forced_src' catalog the source that are in the deepCoad catalogs."""
        deepcoadd = [cat for cat in self.catalogs if 'deepCoadd' in cat]
        if len(deepcoadd):
            if 'forced_src' in self.catalogs:
                print(
                    colored("\nINFO: Matching 'forced_src' and 'deepCoadd' catalogs", "green"))
                print("  - %i sources in the forced-src catalog before selection" %
                      len(self.catalogs['forced_src']))
                coaddid = 'id' if 'id' in self.catalogs[deepcoadd[0]].keys(
                ) else 'objectId'
                filt = numpy.where(numpy.in1d(self.catalogs['forced_src']['objectId'],
                                        self.catalogs[deepcoadd[0]][coaddid]))[0]
                self.catalogs['forced_src'] = self.catalogs['forced_src'][filt]
                print("  - %i sources in the forced-src catalog after selection" %
                      len(self.catalogs['forced_src']))
            else:
                print(colored("\nWARNING: forced_src catalogs not loaded. No match possible.",
                              "yellow"))
        else:
            print(colored("\nWARNING: No deepCoadd* catalog loaded. No match possible.",
                          "yellow"))

    def _add_new_columns(self, catalog=None):
        """Compute magns for all fluxes of a given table. Add the corresponding new columns.

        Compute the x/y position in pixel for all sources. Add new columns to the table.
        """
        print(colored("\nINFO: Adding magnitude and coordinates columns", "green"))
        catalogs = [catalog] if catalog is not None else list(self.catalogs)
        for catalog in catalogs:
            # skip wcs key
            if catalog == 'wcs':
                continue
            print("  - for", catalog)
            columns = []
            # Add magnitudes
            if self.from_butler['getmag'] is not None:
                kfluxes = [
                    k for k in self.catalogs[catalog].columns if k.endswith('_flux')]
                ksigmas = [k + 'Sigma' for k in kfluxes]
                print("    -> getting magnitudes")

                for kflux, ksigma in zip(kfluxes, ksigmas):
                    if kflux.replace('_flux', '_mag') in self.catalogs[catalog].keys():
                        continue

                    if ksigma in self.catalogs[catalog].keys():
                        mag, dmag = self.from_butler['getmag'](numpy.array(self.catalogs[catalog][kflux],
                                                                        dtype='float'),
                                                               numpy.array(self.catalogs[catalog][ksigma],
                                                                        dtype='float'))
                        columns.append(Column(name=kflux.replace('_flux', '_mag'),
                                          data=mag, description='Magnitude', unit='mag'))
                   
                        columns.append(Column(name=ksigma.replace('_fluxSigma', '_magSigma'),
                                          data=dmag, description='Magnitude error', unit='mag'))
            if 'x_Src' in self.catalogs[catalog].keys():
                return
            ra = Quantity(self.catalogs[catalog]["coord_ra"].tolist(), 'rad')
            dec = Quantity(self.catalogs[catalog]["coord_dec"].tolist(), 'rad')
            # Get the x / y position in pixel
            if self.from_butler['wcs'] is not None:
                print("    -> getting pixel coordinates")
                xsrc, ysrc = SkyCoord(ra, dec).to_pixel(
                    self.from_butler['wcs'])
                columns.append(Column(name='x_Src', data=xsrc,
                                      description='x coordinate', unit='pixel'))
                columns.append(Column(name='y_Src', data=ysrc,
                                      description='y coordinate', unit='pixel'))
            else:
                print(colored("\nWARNING: no WCS found for this dataset", "yellow"))

            # Get coordinates in degree
            print("    -> getting degree coordinates")
            columns.append(Column(name='coord_ra_deg', data=Angle(ra).degree,
                                  description='RA coordinate', unit='degree'))
            columns.append(Column(name='coord_dec_deg', data=Angle(dec).degree,
                                  description='DEC coordinate', unit='degree'))

            # Adding all new columns
            print("    -> adding all the new columns")
            self.catalogs[catalog].add_columns(columns)
            # Clean memory before going further
            # gc.collect()

    def _load_calexp(self, calcat='deepCoadd_calexp', **kwargs):
        """Load the deepCoadd_calexp info in order to get the WCS and the magnitudes."""
        print(colored("\nINFO: Loading the %s info" % calcat, 'green'))
        self._load_dataids(calcat, **kwargs)
        print("INFO: Getting the %s catalog for one dataId" % calcat)
        calexp = self._load_catalog_dataid(
            calcat, self.dataids[calcat][0], table=False)
        print("INFO: Getting the magnitude function")
        calib = calexp.getCalib()
        calib.setThrowOnNegativeFlux(False)
        self.from_butler['getmag'] = calib.getMagnitude
        print("INFO: Getting the wcs function")
        wcs = calexp.getWcs().getFitsMetadata().toDict()
        self.from_butler['wcs'] = WCS(wcs)
        self.catalogs['wcs'] = Table({k: [wcs[k]] for k in wcs})

    def load_catalogs(self, catalogs, **kwargs):
        """Load a list of catalogs.

        :param str/list catalogs: A catalog name, or a list of catalogs (see below)
        :param dict keys: A dictionnary of keys to load for each catalog

        Available kwargs are:

        :param bool update: Set to True if you want to update an already loaded catalog
        :param bool show: Set to True to get all available keys of a (list of) catalog(s)
        :param bool matchid: Will only keep objects which are in the deepCoad catalogs (to be used
                             when loading the forced_src and deepCoadd catalogs)

        Examples of catalogs that you can load:

         - 'deepCoadd_ref',
         - 'deepCoadd_meas',
         - 'deepCoadd_forced_src',
         - 'deepCoadd_calexp',
         - 'forced_src'
         - 'src'
        """
        if 'show' in kwargs:
            self.show_keys(catalogs)
            return
        keys = {} if 'keys' not in kwargs else kwargs['keys']
        catalogs = [catalogs] if isinstance(catalogs, str) else catalogs
        if any(["deepCoadd" in cat for cat in catalogs]):
            self._load_calexp(**kwargs)
        else:
            self._load_calexp(calcat='calexp', **kwargs)
        for catalog in sorted(catalogs):
            if catalog in self.catalogs and 'update' not in kwargs:
                print(colored("\nWARNING: %s is already loaded. Use 'update' to reload it." %
                              catalog, "yellow"))
                continue
            if 'calexp' in catalog:
                print(colored("\nWARNING: Skipping %s. Not a regular catalog (no schema).\n" %
                              catalog, "yellow"))
                continue
            print(colored("\nINFO: Loading the %s catalog" % catalog, 'green'))
            self.keys[catalog] = keys.get(catalog, "*")
            self._load_catalog(catalog, **kwargs)
        self._match_deepcoadd_catalogs()
        if 'output_name' in kwargs and self.from_butler['wcs'] is not None:
            self.save_catalogs(kwargs['output_name'],
                               'wcs', kwargs.get('overwrite', False))
        print(colored("\nINFO: Done loading the data.", "green"))

    def show_keys(self, catalogs=None):
        """Show all the available keys."""
        if catalogs is None:
            catalogs = [k for k in self.catalogs.keys() if k != 'wcs']
        catalogs = [catalogs] if isinstance(catalogs, str) else catalogs
        if len(catalogs) == 0:
            print(colored("\nWARNING: No catalog loaded nor given.", "yellow"))
            return
        for cat in catalogs:
            if cat not in self.dataids:
                print(colored("\nINFO: Get the available data IDs", "green"))
                self._load_dataids(cat)
            print(
                colored("\nINFO: Available list of keys for the %s catalog" % cat, "green"))
            table = cutils.get_astropy_table(self.butler.get(cat, dataId=self.dataids[cat][0],
                                                             flags=afwtable.SOURCE_IO_NO_FOOTPRINTS),
            keys="*", get_info=True)
            ktable = Table(numpy.transpose([[k, table[k].description, table[k].unit]
                                         for k in sorted(table.keys())]).tolist(),
                           names=["Keys", "Description", "Units"])
            print("  -> %i keys available for %s" % (len(ktable), cat))
            print("  -> All saved in %s_keys.txt" % cat)
            ktable.write("%s_keys.txt" % cat, format='ascii')

    def save_catalogs(self, output_name, catalog=None, overwrite=False, delete_catalog=False):
        """Save the catalogs into an hdf5 file."""
        # Clean memory before saving
        # gc.collect()
        if not output_name.endswith('.hdf5'):
            output_name += '.hdf5'
        print(colored("\nINFO: Saving the catalogs in %s" % output_name, "green"))
        catalogs = [catalog] if catalog is not None else self.catalogs
        for cat in catalogs:
            print("  - saving", cat)
            for k in self.catalogs[cat].keys():
                if isinstance(self.catalogs[cat][k][0], str):
                    self.catalogs[cat].replace_column(
                        k, Column(self.catalogs[cat][k].astype('bytes')))
            if not self.append:
                self.catalogs[cat].write(output_name, path=cat, compression=True,
                                         serialize_meta=True, overwrite=overwrite)
            else:
                self.catalogs[cat].write(output_name, path=cat, compression=True,
                                         serialize_meta=True, append=True)
            if delete_catalog and cat is not 'wcs':
                oid = self.catalogs[cat]['id' if 'id' in self.catalogs[cat].keys()
                                         else 'objectId'].copy()
                self.catalogs.pop(cat)
                self.catalogs[cat] = Table([oid]).copy()
            self.append = True
        print("INFO: Saving done.")
        # Clean memory before loading a new catalog
        # gc.collect()