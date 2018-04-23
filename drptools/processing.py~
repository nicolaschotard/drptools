import os
import re
import glob
import numpy
import lsst.daf.persistence as dafPersist


DRPPATH = os.getenv('DRPPATH')


class Butler(object):
    
    def __init__(self, drp_path=None):
        
        # Make sure we have data to load
        if drp_path is None and DRPPATH is None:
            raise IOError("You must give a path a DRP output directory.")

        # Load the bulter
        self.drp_path = drp_path if drp_path is not None else DRPPATH
        self.butler = dafPersist.Butler(self.drp_path)
        self.mapper = self.butler._getDefaultMapper()
        self.repoData = self.butler._repos.outputs()[0]
        
        # Load some basic info on the current DRP
        self.repo_input = self._get_repo("input")    # repo containing the rwe data after ingestion
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
                return os.path.realpath(parentRepoData.cfgRoot)                # input path -> ../input in this case
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
    
#    def _get_visit_datasetTypes(self):
#        vds = []
#        for dataset in but.datasetTypes:
#        try:
#            if but.butler.datasetExists(dataset, dataId=but.dataIds['raw'][0]):
#                vds.append(dataset)
#        except:
#            continue

    def get_catIdKeys(self, datasetType):
        """Get the list of ID keys for a given catalog."""
        if datasetType not in self.datasetTypes:
            raise IOError("%s is not a valid datasetType. Check self.datasetTypes to get the valid list." % datasetType) 
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
            return [dict(zip(sorted(keys.keys()), list(v) if not isinstance(v, list) else v)) for v in metadata]
        else:
            if datasetType not in self.repoData.repo._mapper.datasets:
                return []
            template = self.repoData.repo._mapper.datasets[datasetType]._template
            path = os.path.join(self.repoData.cfgRoot, os.path.dirname(template))
            basepath = "/".join([p for p in path.split('/') if not p.startswith('%')]) + "/"
            keys = self.butler.getKeys(datasetType)
            gpath = "/".join([p if not p.startswith('%') else '*' for p in path.split('/')])
            paths = [p for p in glob.glob(gpath) if 'merged' not in p]
            return [{k: keys[k](v) for k, v in zip(keys, p.split(basepath)[1].split('/'))} for p in paths]

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
    
        
