// set Astroem object with important properties and wrapped calls
Astroem = {
      initwcs:  Module.cwrap('initwcs', 'number', ['string', 'number']),
      freewcs:  Module.cwrap('freewcs', 'number', ['number']),
      listhdu:  Module.cwrap('listhdu', 'string', ['string', 'string']),
      wcsinfo:  Module.cwrap('wcsinfo', 'string', ['number']),
      wcssys:  Module.cwrap('wcssys', 'string', ['number', 'string']),
      wcsunits:  Module.cwrap('wcsunits', 'string', ['number', 'string']),
      pix2wcs: Module.cwrap('pix2wcsstr', 'string', ['number', 'number', 'number']),
      wcs2pix: Module.cwrap('wcs2pixstr', 'string', ['number', 'number', 'number']),
      reg2wcs: Module.cwrap('reg2wcsstr', 'string', ['number', 'string']),
      saostrtod: Module.cwrap('saostrtod', 'number', ['string']),
      saodtype:  Module.cwrap('saodtype', 'number'),
      arrfile: Module["arrfile"],
      vfile: Module["vfile"],
      vsize: Module["vsize"],
      vunlink: Module["vunlink"],
      vheap: Module["HEAPU8"],
      vmalloc: Module["_malloc"],
      vmemcpy: Module["writeArrayToMemory"],
      vfree: Module["_free"],
      gzopen: Module.cwrap('gzopen', 'number', ['string', 'string']),
      gzread: Module.cwrap('gzread', 'number', ['number', 'number', 'number']),
      gzwrite: Module.cwrap('gzwrite', 'number', ['number', 'number', 'number']),
      gzclose: Module.cwrap('gzclose', 'number', ['number']),
      gzseek: Module.cwrap('gzseek', 'number', ['number', 'number', 'number']),
      compress: Module['gzcompress'],
      decompress: Module['gzdecompress'],
      handleFITSFile: Module["handleFITSFile"],
      cleanupFITSFile: Module["cleanupFITSFile"],
      getFITSImage: Module["getFITSImage"],
      maxFITSMemory: Module["maxFITSMemory"],
      zscale: Module.cwrap('zscale', 'string', ['number', 'number', 'number', 'number', 'number', 'number', 'number']),
      tanhdr: Module.cwrap('tanhdr', 'string', ['string', 'string', 'string']),
      reproject: Module.cwrap('reproject', 'string', ['string', 'string', 'string', 'string']),
      vls: Module.cwrap('vls', 'int', ['string']),
      vcat: Module.cwrap('vcat', 'string', ['string', 'number']),
      options: Module["options"]
  };

// signal that astroem is ready
if( window.jQuery ){
    if( Module["astroemReady"] ){
	$(document).trigger("astroem:ready", {status: "OK"});
    } else {
	Module["astroemReady"] = true;
    }
}
