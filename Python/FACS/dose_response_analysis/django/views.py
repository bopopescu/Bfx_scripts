# Copyright 2013 XOMA (US) LLC

'''This module defines views.'''

__authors__ = "Mark Evans",


import urllib,os
import pdb
import tempfile, StringIO, zipfile 
from time import localtime,strftime
#from django.contrib.auth.decorators import login_required
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render
from django.core.urlresolvers import reverse


from django import forms
from ri.settings import MEDIA_ROOT
import ri.tools.methods as to_m



#@login_required
def index(request):
    '''Form for uploading FACS MFI data'''
    if request.method == 'POST':
        form = forms.Form(request.POST, request.FILES)
        if form.is_valid():
            filenames = {}
            for upfile in request.FILES.getlist('filename'):
                print "Uploaded ",upfile.name
                # Make sure the facstool subdirectory has been created on the system already
                with tempfile.NamedTemporaryFile(prefix='rawdata_', dir=MEDIA_ROOT+'/facstool/', delete=False) as tmp:
                    content = upfile.read()
                    tmp.write(content)
                    print "temp name = ",tmp.name
                    filenames[tmp.name] = upfile.name
                    tmp.close()
            
            antigens,ag_names = to_m.getAgList(filenames.keys(),filenames)
            request.session['files'] = filenames
            request.session['antigens'] = antigens
            request.session['ag_names'] = ag_names
            request.session.set_expiry(0) # Deletes session data when browser is closed
            url = reverse("tools:defineParameters")

            return HttpResponseRedirect(url)
    else:
        form = forms.Form()
        
    context = { "form": form,
                "form_title": 'Upload Data Files',
                "submit_name": 'Submit',
                "base_template": 'tools/base.html',
                "one_column": True,
                "helptext": "Select one or more data files to upload"}
    
    return render(request, 'tools/data_upload_form.html', context)


def defineParameters(request):

    if request.method =='POST':
        form = forms.Form(request.POST)
        if form.is_valid():
            # Do stuff
            print "Request.POST keys: ",request.POST.keys()
            print "Session keys:",request.session.keys()
            print "Session: ",request.session
            print "Files: ",request.session['files']
            print "Layout: ",request.POST['layout']
            print "Clones: ",request.POST['clones']
            for ag in request.session['ag_names']:
                print "ag renamed: ",ag," -> ",request.POST[ag]
            for m in request.session['ag_names']:
                print "paired ",m," and ",request.POST['paired_'+m]
            for fn in request.session['files'].values():
                print "file= ",fn," condition= ",request.POST[fn+"_cond"]
            ###############
            antigens = request.session['antigens']
            ag_names = request.session['ag_names']
            conc = request.POST['conc']
            log = request.POST['log']
            units = request.POST['units']
            fit = request.POST['fit']
            mod_antigens ={}
            for ag in ag_names:
                if request.POST[ag] == '' or request.POST[ag] == ag: mod_antigens[ag]=ag
                else: mod_antigens[ag] = request.POST[ag]
            files = request.session['files']
            clonelist = request.POST['clones']
            layout = request.POST['layout']
            orig_filenames = antigens.keys()
            conditions = {}
            for ofn in orig_filenames:
                if request.POST[ofn+"_cond"] != '': 
                    antigens[request.POST[ofn+"_cond"]] = antigens[ofn]
                    del(antigens[ofn])
                    for x,y in list(files.iteritems()):
                        if y == ofn:
                            conditions[x] = request.POST[ofn+"_cond"]
            paired_ag = []
            for ag in ag_names:
                if request.POST.has_key("paired_"+ag) and request.POST["paired_"+ag] != '':
                    paired_ag.append((mod_antigens[ag],mod_antigens[request.POST["paired_"+ag]]))

            #############

            #############
            input_file,graphfile,jpg_pathnames = to_m.doDRCAnalysis(antigens,ag_names,files,clonelist,layout,conditions,paired_ag,mod_antigens,conc,units,log,fit)
            
            INFILE = StringIO.StringIO()
            f1 = open(input_file,'r')
            INFILE.write(f1.read())
            f1.close()
            os.remove(input_file)
            PDF = StringIO.StringIO()
            f2 = open(graphfile,'rb')
            PDF.write(f2.read())
            f2.close()
            os.remove(graphfile)
            zip_filename,ZIP = makeZipFile(jpg_pathnames) 
            #os.remove(zip_filename)
            jpg_names = []
            for img in jpg_pathnames:
                _ , img_name=os.path.split(img)
    #            i = open(img,'rb')
     #           IMG = StringIO.StringIO()
      #          IMG.write(i.read())
       #         i.close()
        #        request.session[img_name] = IMG
                jpg_names.append(img_name)
                #os.remove(img)


            request.session['infile'] = INFILE
            request.session['pdffile'] = PDF
            request.session['zipfile'] = ZIP
            context = {'jpg_names':jpg_names}
            return render(request,'tools/show_results.html',context)
            
    else:
        form = forms.Form()

    context = {'form': form,}

    return render(request,'tools/set_parameters.html',context)



def makeZipFile(jpg_names):
    zip_subdir = "drc_graphs"
    zip_filename = "%s.zip" % zip_subdir

    # Open StringIO to grab in-memory ZIP contents
    s = StringIO.StringIO()

    # The zip compressor
    zf = zipfile.ZipFile(s, "w")

    for fpath in jpg_names:
        # Calculate path for file in zip
        fdir, fname = os.path.split(fpath)
        zip_path = os.path.join(zip_subdir, fname)

        # Add file, at correct path
        zf.write(fpath, zip_path)

    # Must close zip for all contents to be written
    zf.close()

    return zip_filename,s


def returnImage(request):
    path = MEDIA_ROOT+"/facstool/"
   # fn = kwargs['fn']
    _, fn = os.path.split(request.META['PATH_INFO'])
    f = open(str(path)+fn,'rb')
    img = f.read()
    f.close()
    os.remove(str(path)+fn)
    response = HttpResponse(img,mimetype='image/jpeg')
    return response


def downloadcsv(request):
    d = strftime("%d%b%Y_%H%M_",localtime())
    tmp = request.session['infile']
    tmp.seek(0)
    response = HttpResponse(tmp.getvalue(),mimetype='text/csv')
    response['Content-Disposition'] = 'attachment;filename='+d+'dose_response_input_data.csv'
    return response


def downloadpdf(request):
    d = strftime("%d%b%Y_%H%M_",localtime())
    tmp = request.session['pdffile']
    tmp.seek(0)
    response = HttpResponse(tmp.getvalue(),mimetype='application/pdf')
    response['Content-Disposition'] = 'attachment; filename='+d+'DRC_graphs.pdf'
    return response


def downloadZipFile(request):
    d = strftime("%d%b%Y_%H%M_",localtime())
    tmp = request.session['zipfile']
    tmp.seek(0)
    # Grab ZIP file from in-memory, make response with correct MIME-type
    response = HttpResponse(tmp.getvalue(), mimetype = "application/x-zip-compressed")
    # ..and correct content-disposition
    response['Content-Disposition'] = ('attachment; filename='+d+'_drc_graphs.zip')
    return response








