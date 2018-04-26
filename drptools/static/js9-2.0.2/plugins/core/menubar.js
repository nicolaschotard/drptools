// ---------------------------------------------------------------------
// JS9 menubar to manage menubar and its menus
// ---------------------------------------------------------------------

/*jslint bitwise: true, plusplus: true, sloppy: true, vars: true, white: true, browser: true, devel: true, continue: true, unparam: true, regexp: true */
/*global $, JS9, sprintf */

// create our namespace, and specify some meta-information and params
JS9.Menubar = {};
JS9.Menubar.CLASS = "JS9";
JS9.Menubar.NAME = "Menubar";
JS9.Menubar.WIDTH = JS9.WIDTH || 512;	// width of Menubar
JS9.Menubar.HEIGHT = "auto";		// height of Menubar

// menu buttons in the menubar
// NB: names must match actual menus, menu labels are arbitrary
JS9.Menubar.buttonOptsArr = [{name: "file", label: "File"},
			     {name: "view", label: "View"},
			     {name: "zoom", label: "Zoom"},
			     {name: "scale", label: "Scale"},
			     {name: "color", label: "Color"},
			     {name: "region", label: "Regions"},
			     {name: "wcs", label: "WCS"},
			     {name: "analysis", label: "Analysis"},
			     {name: "help", label: "Help"}];

// backwards compatibility, pre-1.10
if( JS9.menuButtonOptsArr ){
    JS9.Menubar.buttonOptsArr = JS9.menuButtonOptsArr;
}

// get displays associated with this menubar, taking supermenus into account
JS9.Menubar.getDisplays = function(mode){
    var i, s, disp;
    var arr = [];
    mode = mode || "any";
    if( this.id.search(JS9.SUPERMENU) >= 0 ){
	if( mode !== "all" && this.selectedDisplay ){
	    return [this.selectedDisplay];
	}
	s = this.divjq.data("displays").split(",");
	if( s[0] === "*" ){
	    for(i=0; i<JS9.displays.length; i++){
		arr.push(JS9.displays[i]);
	    }
	} else {
	    for(i=0; i<s.length; i++){
		disp = JS9.lookupDisplay(s[i]);
		if( disp ){
		    arr.push(disp);
		}
	    }
	}
    }
    if( !arr.length ){
	arr.push(this.display);
    }
    return arr;
};

// this callback happens when a click is registered on a display
// we then go through the supermenus, and if one of them contains this display,
// we set its selectedDisplay value so that use of that supermenu is then aimed
// only at the selected display
// also used to unset previously set selectedDisplay
//
// called by JS9.mouseupCB with no context, passing current image object
JS9.Menubar.onclick = function(disp){
    var i, arr, supermenu;
    if( (typeof disp === "string") && (disp !== "all") ){
	disp = JS9.lookupDisplay(disp);
    }
    for(i=0; i<JS9.supermenus.length; i++){
	supermenu = JS9.supermenus[i];
	arr = JS9.Menubar.getDisplays.call(supermenu, "all");
	if( ($.inArray(disp, arr) >= 0) || (disp === "all") ){
	    if( JS9.bugs.webkit_resize ){
		$(".JS9").find(".JS9Image").removeClass("JS9Highlight");
	    } else {
		$(".JS9").removeClass("JS9Highlight");
	    }
	    if( (disp === supermenu.selectedDisplay) || (disp === "all") ){
		// unselect
		supermenu.selectedDisplay = null;
	    } else {
		// select
		supermenu.selectedDisplay = disp;
		if( JS9.bugs.webkit_resize ){
		    $(disp.divjq).find(".JS9Image").addClass("JS9Highlight");
		} else {
		    $(disp.divjq).addClass("JS9Highlight");
		}
	    }
	}
    }
};

JS9.Menubar.init = function(width, height){
    var ii, jj, ss, tt, menu, html;
    var that = this;
    var issuper = this.id.search(JS9.SUPERMENU) >= 0;
    // save object in super array, if necessary
    if( issuper ){
	JS9.supermenus.push(this);
    }
    // set width and height on div
    this.width = this.divjq.attr("data-width");
    if( !this.width  ){
        this.width = width || JS9.Menubar.WIDTH;
        if( issuper ){
            this.width += 10;
        }
    }
    this.divjq.css("width", this.width);
    this.width = parseInt(this.divjq.css("width"), 10);
    this.height = this.divjq.attr("data-height");
    this.buttonClass = this.divjq.attr("data-buttonClass") || "JS9Button" ;
    this.containerClass = "JS9MenubarContainer";
    // special handling of some known button classes
    if( this.buttonClass.match(/-flat/) ){
	this.containerClass += "-flat";
    } else if( this.buttonClass.match(/-border/) ){
	this.containerClass += "-border";
    }
    this.backgroundColor = this.divjq.attr("data-backgroundColor");
    if( !this.height  ){
	this.height = height || JS9.MENUHEIGHT;
    }
    this.divjq.css("height", this.height);
    this.height = parseInt(this.divjq.css("height"), 10);
    // generate html for this menubar
    html = "<span id='JS9Menus_@@ID@@'>";
    if( issuper ){
       html += "<button type='button' id='superMenu@@ID@@'class='"+ this.buttonClass+"'>Super</button>";
    }
    for(jj=0; jj<JS9.globalOpts.menuBar.length; jj++){
	menu = JS9.globalOpts.menuBar[jj];
	for(ii=0; ii<JS9.Menubar.buttonOptsArr.length; ii++){
	    ss = JS9.Menubar.buttonOptsArr[ii].name;
	    if( menu === ss ){
		tt = JS9.Menubar.buttonOptsArr[ii].label;
		// no help available for all-in-one configuration
		if( JS9.allinone && (ss === "help") ){
		    break;
		}
		if( ss[0] === "#" ){
		    ss = ss.slice(1);
		    html += "<button type='button' id='"+ss+"Menu@@ID@@' class='"+ this.buttonClass +"' disabled='disabled'>"+tt+" </button>";
		} else {
		    html += "<button type='button' id='"+ss+"Menu@@ID@@' class='"+ this.buttonClass +"'>"+tt+"</button>";
		}
		break;
	    }
	}
    }
    html += "<button type='button' id='hiddenRegionMenu@@ID@@'class='JS9Button' style='display:none'>R</button>";
    html += "<button type='button' id='hiddenAnchorMenu@@ID@@'class='JS9Button' style='display:none'>R</button>";
    html += "</span>";
    // set the display for this menubar
    this.display = JS9.lookupDisplay(this.id);
    // link back the menubar in the display
    this.display.menubar = this;
    // define menubar
    this.html = html.replace(/@@ID@@/g,this.id);
    // add container to the high-level div
    this.menuConjq = $("<div>")
	.addClass(this.containerClass)
	.attr("width", this.width)
	.attr("height", this.height)
	.html(this.html)
	.appendTo(this.divjq);
    // menubar background color
    if( this.backgroundColor ){
	this.menuConjq.css("background", this.backgroundColor);
    }
    $(function(){
	function onhide() {
	    var tdisp = that.display;
	    if( JS9.bugs.hide_menu && tdisp.image ){
		tdisp.image.displayImage("rgb");
	    }
	}
	if( issuper ){
	// supermenu: make button open the contextMenu
	$("#superMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#superMenu" + that.id).contextMenu();
	});
	$.contextMenu({
            selector: "#superMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, name, tdisp;
		var n = 0;
		var items = {};
		var arr = JS9.Menubar.getDisplays.call(that, "all");
		items.supertitle = {name: "selected display:", disabled: true};
		for(i=0; i<arr.length; i++){
		    tdisp = arr[i];
		    name = tdisp.id;
		    items[name] = {name: name};
		    if( that.selectedDisplay === tdisp ){
			items[name].icon = "sun";
			n++;
		    }
		}
		name = "all displays";
		items.all = {name: name};
		if( !n ){
		    items.all.icon = "sun";
		}
		items["sep" + n++] = "------";
		return{
		    callback: function(key){
			JS9.Menubar.onclick.call(that, key);
		    },
		    items: items
		};
	    }
	});
	}
	// file: make button open the contextMenu
	$("#fileMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#fileMenu" + that.id).contextMenu();
	});
	$.contextMenu({
            selector: "#fileMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, im, name, imlen, s1;
		var n = 0;
		var items = {};
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var tim = tdisp.image;
		items.filetitle = {name: "Images:", disabled: true};
		imlen = JS9.images.length;
		for(i=0; i<imlen; i++){
		    im = JS9.images[i];
		    if( im.display === tdisp ){
			name = im.id;
			if( JS9.globalOpts.rgb.active ){
			    if( im === JS9.globalOpts.rgb.rim){
				name += " (red)";
			    }
			    if( im === JS9.globalOpts.rgb.gim){
				name += " (green)";
			    }
			    if( im === JS9.globalOpts.rgb.bim){
				name += " (blue)";
			    }
			}
			items[name] = {name: name};
			if( tdisp.image && (tdisp.image.id === im.id) ){
			    items[name].icon = "sun";
			}
			n++;
		    }
		}
		if( !n ){
		    items.noimg = {
			name: "[no images]",
			events: {keyup: function(){return;}}
		    };
		}
		items["sep" + n++] = "------";
		items.open = {name: "open local file ..."};
		items.archives = {name: " accessing data archives ..."};
		if( !JS9.allinone ){
		    items.archives.disabled = false;
		} else {
		    items.archives.disabled = true;
		}
		items.loadproxy = {name: "open link via proxy ..."};
		if( !JS9.allinone			 &&
		    JS9.globalOpts.helperType !== "none" &&
		    JS9.globalOpts.workDir      	 &&
		    JS9.globalOpts.loadProxy    	 ){
		    items.loadproxy.disabled = false;
		} else {
		    items.loadproxy.disabled = true;
		}
		items.loadcors = {name: "open link via CORS ..."};
		if( !window.hasOwnProperty("Jupyter") ){
		    items.loadcors.disabled = false;
		} else {
		    items.loadcors.disabled = true;
		}
		items["sep" + n++] = "------";
		items.print = {name: "print ..."};
		items.header = {name: "display FITS header"};
		items.hdus = {name: "display FITS HDUs"};
		if( !tim || !tim.hdus ){
		    items.hdus.disabled = true;
		}
		items.saveas = {
		    name: "save image as ...",
		    items: {saveastitle: {name: "choose output format:",
					  disabled: true},
			    savefits: {name: "FITS"},
			    savejpeg: {name: "JPEG"},
			    savepng: {name: "PNG"}
			   }
		};
		items.moveto = {
		    name: "move this image to ...",
		    items: {movetotitle: {name: "choose display:",
					  disabled: true}}
		};
		items.separate = { name: "separate these images" };
		items.gather = { name: "gather all images here" };
		if( tim ){
		    items.moveto.disabled = false;
		    for(i=0; i<JS9.displays.length; i++){
			if( $("#"+JS9.displays[i].id).length > 0 &&
			    tdisp !== JS9.displays[i]    	     ){
			    s1 = "moveto_" + JS9.displays[i].id;
			    items.moveto.items[s1] = {name: JS9.displays[i].id};
			}
		    }
		    items.moveto.items.moveto_newdisp = { name: "a new display" };
		} else {
		    items.moveto.disabled = true;
		}
		items.free = {name: "free image memory"};
		items.close = {name: "close image"};
		items.closeall = {name: "close all images"};
		items.removeproxy = {name: "remove proxy file from server"};
		if( !tim || !tim.proxyFile ){
		    items.removeproxy.disabled = true;
		}
		items["sep" + n++] = "------";
		items.loadcatalog = {name: "load catalog ..."};
		items.savecatalog = {name: "save active catalog"};
		items["sep" + n++] = "------";
		items.loadsession = {name: "load session ..."};
		items.savesession = {name: "save session"};
		items["sep" + n++] = "------";
		items.lite = {name: "new JS9 light window"};
		items.xnew = {name: "new JS9 separate window"};
		if( window.isElectron ){
		    items.xnew.disabled = true;
		}
		items["sep" + n++] = "------";
		if( window.isElectron && window.electronIPC ){
		    items.electronHelper = {name: "connect to JS9 helper"};
		    if(  JS9.helper.connected ){
			items.electronHelper.disabled = true;
		    }
		}
		items.pageid = {name: "display page id"};
		if( JS9.DEBUG > 2 ){
		  items["sep" + n++] = "------";
		  items.refresh = {name: "debug: refresh ..."};
		}
		return {
                    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var j, s, t, did, kid, unew, uwin;
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			switch(key){
			case "free":
			    if( uim && uim.raw.hdu && uim.raw.hdu.fits ){
				JS9.cleanupFITSFile(uim.raw.hdu.fits, true);
			    }
			    break;
			case "close":
			    if( uim ){
				uim.closeImage();
			    }
			    break;
			case "closeall":
			    if( udisp ){
				// reverse loop because we slice JS9.images
				for(j=JS9.images.length-1; j>=0; j--){
				    uim = JS9.images[j];
				    if( uim.display === udisp ){
					uim.closeImage();
				    }
				}
			    }
			    break;
			case "removeproxy":
			    if( uim ){
				uim.removeProxyFile();
			    }
			    break;
			case "savesession":
			    if( uim ){
				uim.saveSession(uim.id + ".ses");
			    }
			    break;
			case "loadsession":
			    if( udisp ){
				JS9.OpenSessionMenu({display: udisp});
			    }
			    break;
			case "loadcatalog":
			    if( udisp ){
				JS9.OpenCatalogsMenu({display: udisp});
			    }
			    break;
			case "savecatalog":
			    if( uim ){
				uim.saveCatalog();
			    }
			    break;
			case "header":
			    if( uim ){
				if( uim.raw.header ){
				    uim.displayAnalysis("text",
						   JS9.raw2FITS(uim.raw, true),
						   {title: "FITS Header: "+uim.id});
				} else {
				    JS9.error("no FITS header for " + uim.id);
				}
			    }
			    break;
			case "hdus":
			    if( uim ){
				if( uim.hdus ){
				    uim.displayAnalysis("text",
						   JS9.hdus2Str(uim.hdus),
						   {title: "FITS HDUs: "+uim.id,
						    winformat: "width=800px,height=200px,center=1,resize=1,scrolling=1"});
				} else {
				    JS9.error("no FITS header for " + uim.id);
				}
			    }
			    break;
			case "lite":
			    JS9.LoadWindow(null, {clone: udisp.id}, "light");
			    break;
			case "xnew":
			    JS9.LoadWindow(null, null, "new");
			    break;
			case "electronHelper":
			    // Electron.js: send message to main
			    if( window.isElectron && window.electronIPC ){
				try{ window.electronIPC.send("msg",
							     "startHelper"); }
				catch(ignore){}
			    }
			    break;
			case "pageid":
			    s = sprintf("<center><p>pageid: %s</center>",
					JS9.helper.pageid || "none");
			    t = "JS9 page id";
			    // add display to title
			    t += sprintf(JS9.IDFMT, udisp.id);
			    JS9.lightWin("fileid" + JS9.uniqueID(),
					 "inline", s, t,
					 JS9.lightOpts[JS9.LIGHTWIN].lineWin);
			    break;
			case "open":
			    JS9.OpenFileMenu({display: udisp});
			    break;
			case "loadcors":
			    if( JS9.allinone ){
				did = JS9.Image.prototype.displayAnalysis.call(
				      null,
				      "textline",
				      JS9.allinone.loadCorsHTML,
				      {title: "Open a shared CORS link"});
			    } else {
				did = JS9.Image.prototype.displayAnalysis.call(
				      null,
				      "textline",
				      JS9.InstallDir(JS9.globalOpts.corsURL),
				      {title: "Open a shared CORS link"});
			    }
			    // save display id
			    $(did).data("dispid", udisp.id);
			    break;
			case "archives":
			    JS9.DisplayHelp(JS9.InstallDir(JS9.globalOpts.archivesURL));
			    break;
			case "loadproxy":
			    // load param url to run analysis task
			    // param url is relative to js9 install dir
			    did = JS9.Image.prototype.displayAnalysis.call(null,
				     "textline",
				     JS9.InstallDir(JS9.globalOpts.proxyURL),
				     {title: "Open a link via server proxy"});
			    // save info for running the task
			    $(did).data("dispid", udisp.id)
				  .data("aname", "loadproxy");
			    break;
			case "refresh":
			    $('#refreshLocalFile-' + udisp.id).click();
			    break;
			case "savefits":
			    if( uim ){
				s = uim.id.replace(/\.png/i, ".fits")
				          .replace(/\.gz$/i, "")
				          .replace(/\[.*\]/,"");
				uim.saveFITS(s);
			    }
			    break;
			case "savepng":
			    if( uim ){
				s = uim.id.replace(/\.fit[s]?/i, ".png")
				          .replace(/\.gz$/i, "")
				          .replace(/\[.*\]/,"");
				uim.savePNG(s);
			    }
			    break;
			case "savejpeg":
			    if( uim ){
				s = uim.id.replace(/\.fit[s]?/i, ".jpeg")
				          .replace(/\.png$/i, ".jpeg")
				          .replace(/\.gz$/i, "")
				          .replace(/\[.*\]/,"");
				uim.saveJPEG(s);
			    }
			    break;
			case "print":
			    if( uim ){
				uim.print();
			    }
			    break;
			case "separate":
			    if( udisp ){
				udisp.separate();
			    }
			    break;
			case "gather":
			    if( udisp ){
				if( (that.id.search(JS9.SUPERMENU) >= 0) &&
				    !that.selectedDisplay 		 ){
				    JS9.error("gather requires a selected display");
				}
				udisp.gather();
			    }
			    break;
			default:
			    // maybe it's a moveto request
			    if( key.match(/^moveto_/) ){
				unew = key.replace(/^moveto_/,"");
				if( unew === "newdisp" ){
				    uwin = "JS9_light" + JS9.uniqueID();
			            $("#dhtmlwindowholder").arrive("#" + uwin,
                                    {onceOnly: true}, function(){
					uim.moveToDisplay(uwin);
				    });
				    JS9.LoadWindow(null,
                                                   {id: uwin, clone: udisp.id},
                                                   "light");
				} else {
				    uim.moveToDisplay(unew);
				}
				return;
			    }
			    for(j=0; j<JS9.images.length; j++){
				uim = JS9.images[j];
				kid = key.replace(/ *\((red|green|blue)\)/,"");
				if( (udisp.id === uim.display.id) &&
				    (uim.id === kid) ){
				    // display image, 2D graphics, etc.
				    uim.displayImage("all");
				    uim.refreshLayers();
				    udisp.clearMessage();
				    break;
				}
			    }
			    break;
			}
		    });
                    },
		    items: items
		};
            }
	});
	// View: make button open the contextMenu
	$("#viewMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#viewMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#viewMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, plugin, pname, pinst, key;
		var lastxclass="";
		var n = 0;
		var items = {};
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var tim = tdisp.image;
		var editResize = function(disp, obj){
		    var v1, v2, arr;
		    delete tdisp.tmp.editingMenu;
		    if( obj.resize ){
			arr = obj.resize.split(/[\s,\/]+/);
			switch(arr.length){
			case 0:
			    break;
			case 1:
			    if( tim ){
				v1 = tim.wcs2imlen(arr[0]);
				disp.resize(v1, v1);
			    } else if( JS9.isNumber(arr[0]) ){
				v1 = parseInt(arr[0], 10);
				disp.resize(v1, v1);
			    }
			    break;
			default:
			    if( tim && tim.wcs ){
				v1 = tim.wcs2imlen(arr[0]);
				v2 = tim.wcs2imlen(arr[1]);
				disp.resize(v1, v2);
			    } else if( JS9.isNumber(arr[0]) && 
				       JS9.isNumber(arr[1]) ){
				v1 = parseInt(arr[0], 10);
				v2 = parseInt(arr[1], 10);
				disp.resize(v1, v2);
			    }
			    break;
			}
		    }
		};
		var keyResize = function(e){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		    var obj = $.contextMenu.getInputValues(e.data);
		    var keycode = e.which || e.keyCode;
		    var vdisp = val;
		    // make sure display is still valid
		    if( $.inArray(vdisp, JS9.displays) < 0 ){
			return;
		    }
		    switch( keycode ){
		    case 9:
		    case 13:
			editResize(vdisp, obj);
			break;
		    default:
			vdisp.tmp.editingMenu = true;
			break;
		    }
		    });
		};
		// plugins
		items["sep" + n++] = {name: "Plugins:"};
		items["sep" + (n-1)].disabled = true;
		for(i=0; i<JS9.plugins.length; i++){
		    plugin = JS9.plugins[i];
		    pname = plugin.name;
		    if( plugin.opts.menuItem && (plugin.opts.menu === "view") ){
			pinst = tdisp.pluginInstances[pname];
			if( !pinst || pinst.winHandle ){
			    if( plugin.xclass !== lastxclass ){
				// items["sep" + n] = "------";
				n = n + 1;
			    }
			    lastxclass = plugin.xclass;
			    items[pname] = {name: plugin.opts.menuItem};
			    if( pinst && (pinst.status === "active") ){
				items[pname].icon = "sun";
			    }
			}
		    }
		}
		items["sep" + n++] = "------";
		items.valpos = {name: "display value/position"};
		// disable if we don't have info plugin
		if( !JS9.hasOwnProperty("Info") ){
		    items.valpos.disabled = true;
		} else if( tdisp.image && tdisp.image.params.valpos ){
		    items.valpos.icon = "sun";
		}
		items.inherit = {name: "new image inherits current params"};
		if( tdisp.image && tdisp.image.params.inherit ){
		    items.inherit.icon = "sun";
		}
		items["sep" + n++] = "------";
		items.rawlayer = {
		    name: "raw data layers",
		    items: {}
		};
		if( tim && tim.raws.length > 1 ){
		    for(i=0; i<tim.raws.length; i++){
			key = "rawlayer_" + tim.raws[i].id;
			items.rawlayer.items[key] = {
			    name: tim.raws[i].id
			};
			if( tim.raw === tim.raws[i] ){
			    items.rawlayer.items[key].icon = "sun";
			}
		    }
		    items.rawlayer.items["sep" + n++] = "------";
		    items.rawlayer.items.rawlayer_remove = {name: "remove"};
		} else {
		    items.rawlayer.disabled = true;
		}
		items["sep" + n++] = "------";
		items.resize = {
		    events: {keyup: keyResize},
		    name: "change width/height:",
		    type: "text"
		};
		items.imagesize = {name: "set to image size"};
		items.fullsize = {name: "set size to full window"};
		items.resetsize = {name: "reset to original size"};
		if( !JS9.globalOpts.resize ){
		    items.resize.disabled = true;
		    items.fullsize.disabled = true;
		    items.imagesize.disabled = true;
		    items.resetsize.disabled = true;
		} else if( !tim ){
		    items.imagesize.disabled = true;
		}
		return {
		    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		        var jj, ucat, umode, uplugin, s;
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			switch(key){
			case "valpos":
			    if( uim ){
				uim.params.valpos = !uim.params.valpos;
				if( !uim.params.valpos ){
				    udisp.clearMessage();
				}
			    }
			    break;
			case "inherit":
			    if( uim ){
				uim.params.inherit = !uim.params.inherit;
			    }
			    break;
			case "show":
			case "hide":
			    if( uim ){
				for( ucat in uim.layers ){
				    if( uim.layers.hasOwnProperty(ucat) ){
					if( uim.layers[ucat].dlayer.dtype === "main" ){
					    uim.showShapeLayer(ucat, key);
					    if( key === "show" ){
						uim.refreshLayers();
					    }
					}
				    }
				}
			    }
			    break;
			case "fullsize":
			    udisp.resize("full", {center: true});
			    break;
			case "imagesize":
			    udisp.resize("image");
			    break;
			case "resetsize":
			    udisp.resize("reset");
			    break;
			default:
			    // maybe it's a plugin
			    for(jj=0; jj<JS9.plugins.length; jj++){
				uplugin = JS9.plugins[jj];
				if( uplugin.name === key ){
				    udisp.displayPlugin(uplugin);
				    return;
				}
			    }
			    // maybe it's a shape layer
			    if( uim ){
				for( ucat in uim.layers ){
				    if( uim.layers.hasOwnProperty(ucat) ){
					if( key === ucat ){
					    umode = uim.layers[ucat].show ?
						"hide" : "show";
					    uim.showShapeLayer(ucat, umode);
					    if( umode === "show" ){
						uim.refreshLayers();
					    }
					    return;
					}
				    }
				}
			    }
			    // maybe its a raw data layer
			    if( tim && key.match(/^rawlayer_/) ){
				s = key.replace(/^rawlayer_/, "");
				if( s === "remove" ){
				    for(i=0; i<tim.raws.length; i++){
					if( tim.raw === tim.raws[i] ){
					    tim.rawDataLayer(tim.raw.id,
							     "remove");
					}
				    }
				} else {
				    tim.rawDataLayer(s);
				}
			    }
			    break;
			}
		    });
		    },
		    events: {
			show: function(opt){
			    var udisp = that.display;
			    var obj = {};
			    if( udisp  ){
				obj.resize = sprintf("%d %d",
						     udisp.width, udisp.height);
				$.contextMenu.setInputValues(opt, obj);
				JS9.jupyterFocus(".context-menu-item");
			    }
			},
			hide: function(opt){
			    var obj;
			    var udisp = that.display;
			    if( udisp ){
				// if a key was pressed, do the edit
				if( udisp.tmp.editingMenu ){
				    obj = $.contextMenu.getInputValues(opt);
				    editResize(udisp, obj);
				}
			    }
			}
		    },
		    items: items
		};
	    }
	});
	// Zoom: make button open the contextMenu
	$("#zoomMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#zoomMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#zoomMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, zoom, zoomp, name, name2;
		var n = 0;
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var tim = tdisp.image;
		var editZoom = function(im, obj){
		    delete tdisp.tmp.editingMenu;
		    if( !isNaN(obj.zoom) ){
			im.setZoom(obj.zoom);
		    }
		};
		var keyZoom = function(e){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		    var obj = $.contextMenu.getInputValues(e.data);
		    var keycode = e.which || e.keyCode;
		    var vdisp = val;
		    var vim = vdisp.image;
		    // make sure display is still valid
		    if( $.inArray(vdisp, JS9.displays) < 0 ){
			return;
		    }
		    switch( keycode ){
		    case 9:
		    case 13:
			if( vim ){
			    editZoom(vim, obj);
			}
			break;
		    default:
			vdisp.tmp.editingMenu = true;
			break;
		    }
		    });
		};
		var items = {};
		items.zoomtitle = {name: "Zoom Factors:", disabled: true};
		for(i=JS9.imageOpts.zooms; i>=1; i--){
		    zoom = Math.pow(2,-i);
		    zoomp = Math.pow(2,i);
		    name = sprintf("zoom%s", zoom);
		    name2 = sprintf("zoom 1/%s", zoomp);
		    items[name] = {name: name2};
		    if( tim && (tim.rgb.sect.zoom === zoom) ){
			items[name].icon = "sun";
		    }
		}
		for(i=0; i<=JS9.imageOpts.zooms; i++){
		    zoom = Math.pow(2,i);
		    name = sprintf("zoom%s", zoom);
		    name2 = sprintf("zoom %s", zoom);
		    items[name] = {name: name2};
		    if( tim && (tim.rgb.sect.zoom === zoom) ){
			items[name].icon = "sun";
		    }
		}
		items["sep" + n++] = "------";
		items.zoomiotitle = {name: "Zoom In/Out:", disabled: true};
		items.zoomIn = {name: "zoom in"};
		items.zoomOut = {name: "zoom out"};
		items.zoomToFit = {name: "zoom to fit"};
		items["sep" + n++] = "------";
		items.zoom = {
		    events: {keyup: keyZoom},
		    name: "numeric zoom value:",
		    type: "text"
		};
		items["sep" + n++] = "------";
		items.center = {name: "pan to center"};
		items.reset = {name: "reset zoom/pan"};
		return {
		    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			if( uim ){
			    switch(key){
			    case "zoomIn":
				uim.setZoom("x2");
				break;
			    case "zoomOut":
				uim.setZoom("/2");
				break;
			    case "zoomToFit":
				uim.setZoom("tofit");
				break;
			    case "center":
				uim.setPan();
				break;
			    case "reset":
				uim.setZoom("1");
				uim.setPan();
				break;
			    default:
				// look for a numeric zoom
				if( key.match(/^zoom/) ){
				    uim.setZoom(key.slice(4));
				}
				break;
			    }
			}
		    });
		    },
		    events: {
			show: function(opt){
			    var udisp = that.display;
			    var uim = udisp.image;
			    var obj = {};
			    if( uim  ){
				obj.zoom =
				    String(uim.rgb.sect.zoom);
			    }
			    $.contextMenu.setInputValues(opt, obj);
			    JS9.jupyterFocus(".context-menu-item");
			},
			hide: function(opt){
			    var obj;
			    var udisp = that.display;
			    var uim = udisp.image;
			    if( uim ){
				// if a key was pressed, do the edit
				if( udisp.tmp.editingMenu ){
				    obj = $.contextMenu.getInputValues(opt);
				    editZoom(uim, obj);
				}
			    }
			}
		    },
		    items: items
		};
	    }
	});
	// Scale: make button open the contextMenu
	$("#scaleMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#scaleMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#scaleMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, s1, s2;
		var n = 0;
		var items = {};
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var editScale = function(im, obj){
		    delete tdisp.tmp.editingMenu;
		    if( JS9.isNumber(obj.scalemin) ){
			im.params.scalemin = parseFloat(obj.scalemin);
			im.params.scaleclipping = "user";
		    }
		    if( JS9.isNumber(obj.scalemax) ){
			im.params.scalemax = parseFloat(obj.scalemax);
			im.params.scaleclipping = "user";
		    }
		    im.displayImage("colors");
		};
		var keyScale = function(e){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		    var obj = $.contextMenu.getInputValues(e.data);
		    var keycode = e.which || e.keyCode;
		    var vdisp = val;
		    var vim = vdisp.image;
		    // make sure display is still valid
		    if( $.inArray(vdisp, JS9.displays) < 0 ){
			return;
		    }
		    switch( keycode ){
		    case 9:
		    case 13:
			editScale(vim, obj);
			break;
		    default:
			vdisp.tmp.editingMenu = true;
			break;
		    }
		    });
		};
		items.scaletitle = {name: "Scaling Algorithms:",
				    disabled: true};
		for(i=0; i<JS9.scales.length; i++){
		    s1 = JS9.scales[i];
		    s2 = s1;
		    items[s1] = {name: s2};
		    if( tdisp.image && (tdisp.image.params.scale === s1) ){
			items[s1].icon = "sun";
		    }
		}
		items["sep" + n++] = "------";
		items.scalemin = {
		    events: {keyup: keyScale},
		    name: "low limit for clipping:",
		    type: "text"
		};
		items.scalemax = {
		    events: {keyup: keyScale},
		    name: "high limit for clipping:",
		    type: "text"
		};
		items["sep" + n++] = "------";
		items.dminmax = {
		    name: "set limits to data min/max"
		};
		items.zscale = {
		    name: "set limits to zscale z1/z2"
		};
		items.zmax = {
		    name: "set limits to zscale z1/data max"
		};
		return {
                    callback: function(key, opt){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			if( uim ){
			    switch(key){
			    case "dminmax":
				uim.params.scaleclipping = "dataminmax";
				uim.params.scalemin = uim.raw.dmin;
				uim.params.scalemax = uim.raw.dmax;
				$.contextMenu.setInputValues(opt,
				     {scalemin: String(uim.params.scalemin),
				      scalemax: String(uim.params.scalemax)});
				uim.displayImage("colors");
				return false;
			    case "zscale":
				if( (uim.params.z1 === undefined) || 
				    (uim.params.z2 === undefined) ){
				    uim.zscale(false);
				}
				uim.params.scaleclipping = "zscale";
				uim.params.scalemin = uim.params.z1;
				uim.params.scalemax = uim.params.z2;
				$.contextMenu.setInputValues(opt,
				     {scalemin: String(uim.params.scalemin),
				      scalemax: String(uim.params.scalemax)});
				uim.displayImage("colors");
				return false;
			    case "zmax":
				if( (uim.params.z1 === undefined) ){
				    uim.zscale(false);
				}
				uim.params.scaleclipping = "zmax";
				uim.params.scalemin = uim.params.z1;
				uim.params.scalemax = uim.raw.dmax;
				$.contextMenu.setInputValues(opt,
				     {scalemin: String(uim.params.scalemin),
				      scalemax: String(uim.params.scalemax)});
				uim.displayImage("colors");
				return false;
			    default:
				uim.setScale(key);
				break;
			    }
			}
		    });
		    },
		    events: {
			show: function(opt){
			    var udisp = that.display;
			    var uim = udisp.image;
			    var obj = {};
			    if( uim  ){
				obj.scalemin =
				    String(uim.params.scalemin);
				obj.scalemax =
				    String(uim.params.scalemax);
			    }
			    $.contextMenu.setInputValues(opt, obj);
			    JS9.jupyterFocus(".context-menu-item");
			},
			hide: function(opt){
			    var obj;
			    var udisp = that.display;
			    var uim = udisp.image;
			    if( uim ){
				// if a key was pressed, do the edit
				if( udisp.tmp.editingMenu ){
				    obj = $.contextMenu.getInputValues(opt);
				    editScale(uim, obj);
				}
			    }
			}
		    },
		    items: items
		};
	    }
	});
	// Color: make button open the contextMenu
	$("#colorMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#colorMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#colorMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, s1, s2, arr;
		var n = 0;
		var items = {};
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var editColor = function(im, obj){
		    delete tdisp.tmp.editingMenu;
		    if( obj.contrast && !isNaN(obj.contrast) ){
			im.params.contrast = parseFloat(obj.contrast);
		    }
		    if( obj.bias && !isNaN(obj.bias) ){
			im.params.bias = parseFloat(obj.bias);
		    }
		    if( !isNaN(obj.opacity) ){
			if( obj.opacity !== "" ){
			    im.params.opacity = parseFloat(obj.opacity);
			} else {
			    im.params.opacity = 1.0;
			}
		    }
		    im.displayImage("colors");
		};
		var keyColor = function(e){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		    var obj = $.contextMenu.getInputValues(e.data);
		    var keycode = e.which || e.keyCode;
		    var vdisp = val;
		    var vim = vdisp.image;
		    // make sure display is still valid
		    if( $.inArray(vdisp, JS9.displays) < 0 ){
			return;
		    }
		    switch( keycode ){
		    case 9:
		    case 13:
			editColor(vim, obj);
			break;
		    default:
			vdisp.tmp.editingMenu = true;
			break;
		    }
		    });
		};
		items.cmaptitle = {name: "Colormaps:", disabled: true};
		for(i=0; i<JS9.globalOpts.topColormaps.length; i++){
		    s1 = JS9.globalOpts.topColormaps[i];
		    s2 = s1;
		    items[s1] = {name: s2};
		    if( tdisp.image && (tdisp.image.cmapObj.name === s1) ){
			items[s1].icon = "sun";
		    }
		}
		items.morecmaps = {
		    name: "more ...",
		    items: {morecmapstitle: {name: "Colormaps:",
					  disabled: true}}
		};
		for(i=0; i<JS9.colormaps.length; i++){
		    s1 = JS9.colormaps[i].name;
		    if( JS9.globalOpts.topColormaps.indexOf(s1) === -1 ){
		        s2 = s1;
		        items.morecmaps.items[s1] = {name: s2};
		        if( tdisp.image && (tdisp.image.cmapObj.name === s1) ){
			    items.morecmaps.items[s1].icon = "sun";
		        }
                    }
		}
		items["sep" + n++] = "------";
		items.imfilter = {
		    name: "image filters",
		    items: {}
		};
		arr = JS9.Image.prototype.filterRGBImage.call(null).sort();
		for(i=0; i<arr.length; i++){
		    if( arr[i] === "convolve" ){
			continue;
		    }
		    s1 = "imfilter_" + arr[i];
		    items.imfilter.items[s1] = {
			name: arr[i]
		    };
		}
		items["sep" + n++] = "------";
		items.contrast = {
		    events: {keyup: keyColor},
		    name: "contrast value:",
		    type: "text"
		};
		items.bias = {
		    events: {keyup: keyColor},
		    name: "bias value:",
		    type: "text"
		};
		items.opacity = {
		    events: {keyup: keyColor},
		    name: "opacity value:",
		    type: "text"
		};
		items["sep" + n++] = "------";
		items.reset = {name: "reset contrast/bias"};
		items["sep" + n++] = "------";
		items.loadcmap = {name: "load colormap"};
		items.savecmap = {name: "save colormap"};
		items.invert = {name: "invert colormap"};
		if( tdisp.image && tdisp.image.params.invert ){
		    items.invert.icon = "sun";
		}
		items["sep" + n++] = "------";
		items.rgb = {name: "RGB mode"};
		if( JS9.globalOpts.rgb.active ){
		    items.rgb.icon = "sun";
		}
		return {
		    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			if( uim ){
			    switch(key){
			    case "loadcmap":
				JS9.OpenColormapMenu({display: udisp});
				break;
			    case "savecmap":
				JS9.SaveColormap({display: udisp});
				break;
			    default:
				if( key.match(/^imfilter_/) ){
				    s1 = key.replace(/^imfilter_/,"");
				    uim.filterRGBImage(s1);
				    return;
				}
				uim.setColormap(key);
			    }
			}
		    });
		    },
		    events: {
			show: function(opt){
			    var udisp = that.display;
			    var uim = udisp.image;
			    var obj = {};
			    if( uim  ){
				obj.contrast = String(uim.params.contrast);
				obj.bias = String(uim.params.bias);
				obj.opacity = String(uim.params.opacity);
				obj.sigma = String(uim.params.sigma);
			    }
			    $.contextMenu.setInputValues(opt, obj);
			    JS9.jupyterFocus(".context-menu-item");
			},
			hide: function(opt){
			    var obj;
			    var udisp = that.display;
			    var uim = udisp.image;
			    if( uim ){
				// if a key was pressed, do the edit
				if( udisp.tmp.editingMenu ){
				    obj = $.contextMenu.getInputValues(opt);
				    editColor(uim, obj);
				}
			    }
			}
		    },
		    items: items
		};
	    }
	});
	// Region: make button open the contextMenu
	$("#regionMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#regionMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#regionMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, s1;
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var tim = tdisp.image;
		var items = {
		    "regiontitle": {name: "Regions:", disabled: true},
		    "annulus": {name: "annulus"},
		    "box": {name: "box"},
		    "circle": {name: "circle"},
		    "ellipse": {name: "ellipse"},
		    "line": {name: "line"},
		    "point": {name: "point"},
		    "polygon": {name: "polygon"},
		    "text": {name: "text"},
		    "sep1": "------",
		    "loadRegions" : {name: "load regions"},
		    "saveRegions" : {name: "save regions"},
		    "copyto" : {name: "copy regions to",
				items: { copytotitle: {name: "choose image:",
						       disabled: true}}},
		    "listRegions" : {name: "list regions"},
		    "removeRegions" : {name: "remove regions"},
		    "sep2": "------",
		    "listSelRegions" : {name: "list selected"},
		    "removeSelRegions" : {name: "remove selected"},
		    "sep3": "------",
		    "listonchange" : {name: "list on change"},
		    "xeqonchange" : {name: "xeq on change"}
		};
		if( tim && (JS9.images.length > 1) ){
		    for(i=0; i<JS9.images.length; i++){
			if( tim !== JS9.images[i] ){
			    s1 = "copyto_" + JS9.images[i].id;
			    items.copyto.items[s1] = {name: JS9.images[i].id};
			}
		    }
		    items.copyto.items.copyto_all = { name: "all images" };
		    items.copyto.disabled = false;
		} else {
		    items.copyto.disabled = true;
		}
		// disable if we don't have info plugin
		if( !JS9.hasOwnProperty("Info") ){
		    items.listRegions.disabled = true;
		}
		if( tim && tim.params.listonchange ){
		    items.listonchange.icon = "sun";
		}
		if( tim && tim.params.xeqonchange ){
		    items.xeqonchange.icon = "sun";
		}
		return {
		    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var uid;
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			if( uim ){
			    switch(key){
			    case "removeRegions":
				uim.removeShapes("regions", "all");
				udisp.clearMessage("regions");
				break;
			    case "removeSelRegions":
				uim.removeShapes("regions", "selected");
				udisp.clearMessage("regions");
				break;
			    case "loadRegions":
				JS9.OpenRegionsMenu({display: udisp});
				break;
			    case "saveRegions":
				uim.saveRegions("js9.reg", "all");
				break;
			    case "listRegions":
				uim.listRegions("all", {mode: 2});
				break;
			    case "listSelRegions":
				uim.listRegions("selected", {mode: 2});
				break;
			    case "xeqonchange":
				uim.params.xeqonchange = !uim.params.xeqonchange;
				break;
			    case "listonchange":
				uim.params.listonchange = !uim.params.listonchange;
				break;
			    default:
				// maybe it's a copyto request
				if( key.match(/^copyto_/) ){
				    uid = key.replace(/^copyto_/,"");
				    uim.copyRegions(uid);
				    return;
				}
				// otherwise it's new region
				uim.addShapes("regions", key, {ireg: true});
				break;
			    }
			}
		    });
		    },
		    items: items
		};
	    }
	});
	// WCS: make button open the contextMenu
	$("#wcsMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#wcsMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#wcsMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var i, s1, s2, key, altwcs;
		var n=0, nwcs=0, got=0;
		var items = {};
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var tim = tdisp.image;
		var editRotate = function(im, obj){
		    delete tdisp.tmp.editingMenu;
		    if( JS9.isNumber(obj.rot) ){
			im.rotateData(parseFloat(obj.rot));
		    }
		};
		var keyRotate = function(e){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		    var obj = $.contextMenu.getInputValues(e.data);
		    var keycode = e.which || e.keyCode;
		    var vdisp = val;
		    var vim = vdisp.image;
		    // make sure display is still valid
		    if( $.inArray(vdisp, JS9.displays) < 0 ){
			return;
		    }
		    switch( keycode ){
		    case 9:
		    case 13:
			if( vim ){
			    editRotate(vim, obj);
			}
			break;
		    default:
			vdisp.tmp.editingMenu = true;
			break;
		    }
		    });
		};
		items.wcssystitle = {name: "WCS Systems:", disabled: true};
		for(i=0; i<JS9.wcssyss.length; i++){
		    s1 = JS9.wcssyss[i];
		    s2 = s1;
		    items[s1] = {name: s2};
		    if( tim && (tim.params.wcssys === s1) ){
			items[s1].icon = "sun";
			got++;
		    }
		}
		// if we don't know which wcssys is current, assume "native"
		if( !got ){
		    s1 = "native";
		    items[s1].icon = "sun";
		}
		items["sep" + n++] = "------";
		items.wcsutitle = {name: "WCS Units:", disabled: true};
		for(i=0; i<JS9.wcsunitss.length; i++){
		    s1 = JS9.wcsunitss[i];
		    s2 = s1;
		    items[s1] = {name: s2};
		    if( tim && (tim.params.wcsunits === s1) ){
			items[s1].icon = "sun";
		    }
		}
		items["sep" + n++] = "------";
		items.altwcs = {
		    name: "alternate wcs",
		    items: {altwcstitle: {name: "choose a wcs:",
					  disabled: true}}
		};
		if( !tim || !tim.raw || !tim.raw.altwcs ){
		    items.altwcs.disabled = true;
		} else {
		    altwcs = tim.raw.altwcs;
		    for(key in altwcs ){
			if( altwcs.hasOwnProperty(key) ){
			    s1 = "altwcs_" + key;
			    if( altwcs[key].header.WCSNAME ){
				s2 = altwcs[key].header.WCSNAME + 
				    "    (" + key + ")";
			    } else {
				s2 = key;
			    }
			    items.altwcs.items[s1] = { name: s2 };
			    if( tim.raw.wcs === altwcs[key].wcs ){
				items.altwcs.items[s1].icon = "sun";
			    }
			    nwcs++;
			}
		    }
		    // disable if we only have the default wcs
		    if( nwcs < 2 ){
			items.altwcs.disabled = true;
			items.altwcs.items.notasks = {
			    name: "[none]",
			    disabled: true,
			    events: {keyup: function(){return;}}
			};
		    }
		}
		items["sep" + n++] = "------";
		items.reproject = {
		    name: "wcs reproject",
		    items: {reprojtitle: {name: "using the wcs from:", disabled: true}}
		};
		for(i=0, nwcs=0; i<JS9.images.length; i++){
		    if( JS9.images[i].raw.wcs ){
			if( (tim === JS9.images[i]) &&
			    (that.id.search(JS9.SUPERMENU) < 0) ){
			    continue;
			}
			s1 = "reproject_" + JS9.images[i].id;
			items.reproject.items[s1] = {
			    name: JS9.images[i].id
			};
			nwcs++;
		    }
		}
		if( nwcs === 0 ){
		    items.reproject.items.notasks = {
			name: "[none]",
			disabled: true,
			events: {keyup: function(){return;}}
		    };
		} else {
		    items.reproject.disabled = false;
		    items.reproject.items["sep" + n++] = "------";
		    items.reproject.items.reproject_wcsalign = {
			name: "display wcs-aligned"
		    };
		    if( tim && (tim.params.wcsalign) ){
			items.reproject.items.reproject_wcsalign.icon = "sun";
		    }
		}
		items.reproject.items["sep" + n++] = "------";
		items.reproject.items.rotatetitle = {
		    name: "by rotating the image:", disabled: true
		};
		items.reproject.items.reproject_northup = {
		    name: "so that north is up"
		};
		items.reproject.items.rot = {
		    events: {keyup: keyRotate},
		    name: "using angle in degrees:",
		    type: "text"
		};
		if( !tim || !tim.raw || !tim.raw.header || !tim.raw.wcsinfo ){
		    items.reproject.disabled = true;
		}
		return {
                    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var file, s;
			var rexp = new RegExp(key);
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			if( uim ){
			    // maybe it's an alt wcs request
			    if( key.match(/^altwcs_/) ){
				s = key.replace(/^altwcs_/,"");
				uim.setWCS(s);
				return;
			    }
			    // maybe it's a wcs reprojection request
			    if( key.match(/^reproject_/) ){
				if( key === "reproject_wcsalign" ){
				    uim.params.wcsalign = !uim.params.wcsalign;
				    uim.displayImage("display");
				} else if( key === "reproject_northup" ){
				    uim.rotateData("northisup");
				}  else {
				    file = key.replace(/^reproject_/,"");
				    uim.reprojectData(file);
				}
				return;
			    }
			    // otherwise it's a wcs directive
			    if( JS9.wcssyss.join("@").search(rexp) >=0 ){
				uim.setWCSSys(key);
				uim.updateShapes("regions", "all", "wcs");
			    } else if( JS9.wcsunitss.join("@").search(rexp)>=0){
				uim.setWCSUnits(key);
				uim.updateShapes("regions", "all", "wcs");
			    } else {
				JS9.error("unknown wcs sys/units: " + key);
			    }
			}
		    });
		    },
		    events: {
			show: function(opt){
			    var udisp = that.display;
			    var uim = udisp.image;
			    var obj = {};
			    if( uim ){
				obj.rot = "";
				$.contextMenu.setInputValues(opt, obj);
				JS9.jupyterFocus(".context-menu-item");
			    }
			},
			hide: function(opt){
			    var obj;
			    var udisp = that.display;
			    var uim = udisp.image;
			    if( uim ){
				obj = $.contextMenu.getInputValues(opt);
				// if a key was pressed, do the edit
				if( udisp.tmp.editingMenu ){
				    editRotate(uim, obj);
				}
			    }
			}
		    },
		    items: items
		};
	    }
	});
	// ANALYSIS: make button open the contextMenu
	$("#analysisMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#analysisMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#analysisMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
	        var i, j, s, apackages, atasks;
		var plugin, pinst, pname;
		var parr;
		var parexp = /fitsHeader\(([A-Za-z0-9_]+),(.*)\)/;
		var winexp = /winVar\((.*),(.*)\)/;
		var js9exp = /js9Var\((.*),(.*)\)/;
		var imexp = /imVar\((.*),(.*)\)/;
		var ntask = 0;
		var n = 0;
		// var m = 0;
		var items = {};
		var tdisp = JS9.Menubar.getDisplays.call(that)[0];
		var im = tdisp.image;
		var lastxclass="";
		var seq = function(s1, s2){
		    if( !s1 || !s2 ){
			return false;
		    }
		    return String(s1).toUpperCase() === 
			   String(s2).toUpperCase();
		};
		var editAnalysis = function(im, obj){
		    delete tdisp.tmp.editingMenu;
		    if( !obj.sigma ){
			obj.sigma = "none";
		    }
		    if( obj.sigma === "none" ){
			im.params.sigma = obj.sigma;
			im.rawDataLayer("gaussBlur", "remove");
		    } else {
			im.params.sigma = parseFloat(obj.sigma);
			im.gaussBlurData(im.params.sigma);
		    }
		};
		var keyAnalysis = function(e){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
		    var obj = $.contextMenu.getInputValues(e.data);
		    var keycode = e.which || e.keyCode;
		    var vdisp = val;
		    var vim = vdisp.image;
		    // make sure display is still valid
		    if( $.inArray(vdisp, JS9.displays) < 0 ){
			return;
		    }
		    switch( keycode ){
		    case 9:
		    case 13:
			editAnalysis(vim, obj);
			break;
		    default:
			vdisp.tmp.editingMenu = true;
			break;
		    }
		    });
		};
		for(i=0; i<JS9.plugins.length; i++){
		    plugin = JS9.plugins[i];
		    pname = plugin.name;
		    if( plugin.opts.menuItem &&
			(plugin.opts.menu === "analysis") ){
			pinst = tdisp.pluginInstances[pname];
			if( !pinst || pinst.winHandle ){
			    if( plugin.xclass !== lastxclass ){
				if( n > 0 ){
				    items["sep" + n++] = "------";
				}
				items["sep" + n++] =
				    {name: plugin.xclass + " Plugins:"};
			        items["sep" + (n-1)].disabled = true;
			    }
			    lastxclass = plugin.xclass;
			    items[pname] = {
				name: plugin.opts.menuItem
			    };
			    if( pinst && (pinst.status === "active") ){
				items[pname].icon = "sun";
			    }
			    n++;
			}
		    }
		}
		// no server side analysis for CDN all-in-one configuration
		if( !JS9.allinone ){
		if( n > 0 ){
		    items["sep" + n++] = "------";
		}
	        items.remotetitle = {
		    name: "Server-side Tasks:",
		    disabled: true
	        };
		if( im && im.analysisPackages ){
		    apackages = im.analysisPackages;
		    // m = 0;
		    for(j=0; j<apackages.length; j++){
			atasks = apackages[j];
			for(i=0; i<atasks.length; i++){
			    if( atasks[i].hidden ){
				continue;
			    }
			    if( atasks[i].files.match(/^fits$/) &&
				!im.fitsFile ){
				continue;
			    }
			    if( atasks[i].files.match(/^png$/) &&
				(im.source !== "fits2png") ){
				continue;
			    }
			    if( atasks[i].files.match(/^table$/) ){
				if( im.imtab !== "table" ){
				    continue;
				}
			    }
			    if( atasks[i].files.match(/^image$/) ){
				if( im.imtab !== "image" ){
				    continue;
				}
			    }
			    // header params: fitsHeader(pname,pvalue)
			    parr = atasks[i].files.match(parexp);
			    if( parr ){
				s = im.raw.header[parr[1].toUpperCase()];
				if( !seq(s, parr[2]) ){
				    continue;
				}
			    }
			    // win vars: winVar(name,value)
			    parr = atasks[i].files.match(winexp);
			    if( parr ){
				s = JS9.varByName(parr[1], window);
				if( !seq(s, parr[2]) ){
				    continue;
				}
			    }
			    // js9 vars: js9Var(name,value)
			    parr = atasks[i].files.match(js9exp);
			    if( parr ){
				s = JS9.varByName(parr[1], JS9);
				if( !seq(s, parr[2]) ){
				    continue;
				}
			    }
			    // im vars: imVar(name,value)
			    parr = atasks[i].files.match(imexp);
			    if( parr ){
				s = JS9.varByName(parr[1], im);
				if( !seq(s, parr[2]) ){
				    continue;
				}
			    }
			    s = atasks[i].title;
			    if( atasks[i].purl ){
				s += " ...";
			    }
			    items[atasks[i].name] = {
				name: s
			    };
			    ntask++;
			    // m++;
			}
		    }
		}
		if( !ntask ){
		    items.notasks = {
			name: "[none]",
			disabled: true,
			events: {keyup: function(){return;}}
		    };
		    if( JS9.globalOpts.loadProxy &&
			im && im.raw && im.raw.hdu && im.raw.hdu.vfile ){
			items.upload = {name: "upload FITS to make tasks available"};
			if( !JS9.helper.connected || (JS9.helper.type !== "nodejs" && JS9.helper.type !== "socket.io") ){
			    items.upload.disabled = true;
			}
		    }
		}
		items["sep" + n++] = "------";
		items.sigma = {
		    events: {keyup: keyAnalysis},
		    name: "Gaussian blur sigma:",
		    type: "text"
		};
		items.dpath = {name: "set data path ..."};
		}
		return {
                    callback: function(key){
		    JS9.Menubar.getDisplays.call(that).forEach(function(val){
			var a, did, jj, tplugin;
			var udisp = val;
			var uim = udisp.image;
			// make sure display is still valid
			if( $.inArray(udisp, JS9.displays) < 0 ){
			    return;
			}
			// first look for a plugin -- no image rquired
			for(jj=0; jj<JS9.plugins.length; jj++){
			    tplugin = JS9.plugins[jj];
			    if( tplugin.name === key ){
				udisp.displayPlugin(tplugin);
				return;
			    }
			}
			// the rest need an image loaded
			if( uim ){
			    switch(key){
			    case "dpath":
				// call this once window is loaded
			        $("#dhtmlwindowholder").arrive("#dataPathForm",
                                {onceOnly: true}, function(){
				    $('#dataPath').val(JS9.globalOpts.dataPath);
				});
				did = uim.displayAnalysis("textline",
					 JS9.InstallDir(JS9.analOpts.dpathURL),
					 {title: "Data Path for Drag and Drop"});
				// save display id
				$(did).data("dispid", udisp.id);
				break;
			    case "upload":
				uim.uploadFITSFile();
				break;
			    default:
				// look for analysis routine
				a = uim.lookupAnalysis(key);
				if( a ){
				    // load param url to run analysis task
				    // param url is relative to js9 install dir
				    if( a.purl ){
					did = uim.displayAnalysis("params",
						  JS9.InstallDir(a.purl),
						  {title: a.title+": "+uim.fitsFile,
						  winformat: a.pwin});
					// save info for running the task
					$(did).data("dispid", udisp.id)
				              .data("aname", a.name);
				    } else {
					// else run task directly
					uim.runAnalysis(a.name);
				    }
				}
				return;
			    }
			}
		    });
		    },
		    events: {
			show: function(opt){
			    var udisp = that.display;
			    var uim = udisp.image;
			    var obj = {};
			    if( uim  ){
				obj.sigma = String(uim.params.sigma);
			    }
			    $.contextMenu.setInputValues(opt, obj);
			    JS9.jupyterFocus(".context-menu-item");
			},
			hide: function(opt){
			    var obj;
			    var udisp = that.display;
			    var uim = udisp.image;
			    if( uim ){
				// if a key was pressed, do the edit
				if( udisp.tmp.editingMenu ){
				    obj = $.contextMenu.getInputValues(opt);
				    editAnalysis(uim, obj);
				}
			    }
			}
		    },
		    items: items
		};
	    }
	});
	// HELP: make button open the contextMenu
	$("#helpMenu" + that.id).on("mousedown", function(evt){
            evt.preventDefault();
            $("#helpMenu" + that.id).contextMenu();
	});
	// define contextMenu actions
	$.contextMenu({
            selector: "#helpMenu" + that.id,
	    zIndex: JS9.MENUZINDEX,
	    events: { hide: onhide },
            build: function(){
		var key, val;
		var n = 1;
		var last = "";
		var items = {};
		items.helptitle = {name: "JS9 help:", disabled: true};
		// first, internal js9 pages
		for( key in JS9.helpOpts ){
		    if( JS9.helpOpts.hasOwnProperty(key) ){
			val = JS9.helpOpts[key];
			if( val.heading === "JS9Help" ){
			    last = val.type;
			    items[key] = {name: val.title};
			}
		    }
		}
		items["sep" + n++] = "------";
		items.plugintitle = {name: "JS9 plugins:", disabled: true};
		// second, the JS9 core plugins
		for( key in JS9.helpOpts ){
		    if( JS9.helpOpts.hasOwnProperty(key) ){
			val = JS9.helpOpts[key];
			if( val.heading === "JS9" ){
			    last = val.type;
			    items[key] = {name: val.title};
			}
		    }
		}
		// last, the others
		for( key in JS9.helpOpts ){
		    if( JS9.helpOpts.hasOwnProperty(key) ){
			val = JS9.helpOpts[key];
			if( val.heading === "JS9Help" || 
			    val.heading === "JS9" ){
			    continue;
			}
			if( (last !== "") && (val.type !== last) ){
			    items["sep" + n++] = "------";
			    if( val.heading ){
				items["sep" + n++] =
				    {name: val.heading + " plugins:"};
			        items["sep" + (n-1)].disabled = true;
			    }
			}
			last = val.type;
			items[key] = {name: val.title};
		    }
		}
		items["sep" + n++] = "------";
		items.about = {name: "About JS9"};
		return{
		    callback: function(key){
			switch(key){
			case "about":
			    alert(sprintf("JS9: astronomical image display everywhere\nversion: %s\nprincipals: Eric Mandel (lead), Alexey Vikhlinin (science,management)\ncontact: saord@cfa.harvard.edu\n%s", JS9.VERSION, JS9.COPYRIGHT));
			    break;
			default:
			    JS9.DisplayHelp(key);
			    break;
			}
		    },
		    items: items
		};
	    }
	});
    });
};

JS9.RegisterPlugin("JS9", "Menubar", JS9.Menubar.init,
		   {winDims: [JS9.Menubar.WIDTH, JS9.Menubar.HEIGHT]});

