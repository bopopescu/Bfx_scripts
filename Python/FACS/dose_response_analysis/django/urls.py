# Copyright 2013 XOMA (US) LLC

'''URLs for the Dose Response Curve app.'''

__authors__ = "Mark Evans",


from django.conf.urls import patterns, url
import ri.tools.views as rv

urlpatterns = patterns('ri.tools.views',
    url(r'^$', 'index', name="index"),
    url(r'defineParameters','defineParameters',name='defineParameters'),
    url(r'downloadcsv','downloadcsv',name='downloadcsv'),
    url(r'downloadpdf','downloadpdf',name='downloadpdf'),
    url(r'downloadZipFile','downloadZipFile',name='downloadZipFile'),
   # url(r'returnImage/(?P<fn>.*)',rv.returnImage,name='returnImage'),
    #url(r'returnImage/(?P<fn>.*)',rv.returnImage,name='returnImage'),
    url(r'returnImage/.*','returnImage',name='returnImage'),
   # url(r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT}, name="media")
    )

