#!/usr/bin/env python

### T.Navarro

##
#   This routine reads sparse data, i.e. not binned in a lat/lon/vert/time grid and does plots
#
#   There is nothing gnuplot could not be able to do for now, except command line and ease and all python plots possibilities, 
#   but the goal is to provide a framework to compare with data from simulations on a constant grid.
###

###########################################################################################
###########################################################################################
### What is below relate to running the file as a command line executable (very convenient)
if __name__ == "__main__":
   import sys
   from optparse import OptionParser    ### to be replaced by argparse
   from os import system,path
   import numpy as np
   from sparse import sparse

   parser = OptionParser()

   #############################
   ### Options
   parser.add_option('-f', '--file',    action='store',dest='file',     type="string",  default=None,  help='[NEEDED] filename, comma separated')
   parser.add_option('--ftype',         action='store',dest='ftype',   default=None, help=' Force data file type [\'sav\',\'txt\']')
   parser.add_option('-X', '--xvar',    action='store',dest='xvar',     type="string",  default=None,  help='[NEEDED] x variable: a name, a number, etc ...')
   parser.add_option('-Y', '--yvar',    action='store',dest='yvar',     type="string",  default=None,  help='y variable')
   parser.add_option('-Z', '--zvar',    action='store',dest='zvar',     type="string",  default=None,  help='z variable')
   parser.add_option('-C', '--cvar',    action='store',dest='cvar',     type="string",  default=None,  help='c variable for color')
   parser.add_option('--xmin',          action='store',dest='xmin',     type="float",   default=None,  help='min values for x var, comma separated')
   parser.add_option('--xmax',          action='store',dest='xmax',     type="float",   default=None,  help='max values for x var, comma separated')
   parser.add_option('--ymin',          action='store',dest='ymin',     type="float",   default=None,  help='min values for y var, comma separated')
   parser.add_option('--ymax',          action='store',dest='ymax',     type="float",   default=None,  help='max values for y var, comma separated')
   parser.add_option('--zmin',          action='store',dest='zmin',     type="float",   default=None,  help='min values for z var, comma separated')
   parser.add_option('--zmax',          action='store',dest='zmax',     type="float",   default=None,  help='max values for z var, comma separated')
   parser.add_option('--cmin',          action='store',dest='cmin',     type="string",  default=None,  help='min values for c var, comma separated')
   parser.add_option('--cmax',          action='store',dest='cmax',     type="string",  default=None,  help='max values for c var, comma separated')
   parser.add_option('-w','--with',     action='append',dest='cond',    type="string",  default=None,  help='conditions..')
   parser.add_option('--merge',   action='store_true', dest='merge',                      default=False, help='merge datafiles in a single plot [False]')
   parser.add_option('-c','--color',  action='store',dest='clb',       type="string",  default="def", help='change colormap (also: nobar,onebar)')
   parser.add_option('--trans',        action='store',dest='trans',     type="float",   default=1.,    help='shaded plot transparency, 0 to 1 (=opaque) [1]')
   parser.add_option('-p', '--proj',   action='store',dest='proj',      type="string",  default=None,  help='projection')
   parser.add_option('--blat',         action='store',dest='blat',      type="int",     default=None,  help='reference lat (or bounding lat for stere) [computed]')
   parser.add_option('--blon',         action='store',dest='blon',      type="int",     default=None,  help='reference lon [computed]')
   parser.add_option('-b', '--back',   action='store',dest='back',      type="string",  default=None,  help='background image [None]')
   parser.add_option('-m', '--min',    action='store',dest='vmin',      type="float",   default=None,  help='bounding minimum value [min]')
   parser.add_option('-M', '--max',    action='store',dest='vmax',      type="float",   default=None,  help='bounding maximum value [max]')
   parser.add_option('--div',          action='store',dest='ndiv',      type="int",     default=10,    help='number of divisions in histogram [10]')
   parser.add_option('-S', '--save',   action='store',dest='save',      type="string",  default="gui", help='save mode (gui,png,eps,svg,pdf,txt,html,avi) [gui]')


   #############################
   ### Get options and variables
   (opt,args) = parser.parse_args()

   sparse(oplot=None,file=opt.file,ftype=opt.ftype,xvar=opt.xvar,yvar=opt.yvar,zvar=opt.zvar,cvar=opt.cvar,xmin=opt.xmin,xmax=opt.xmax,ymin=opt.ymin,ymax=opt.ymax,zmin=opt.zmin,zmax=opt.zmax,cmin=opt.cmin,cmax=opt.cmax,cond=opt.cond,merge=opt.merge,clb=opt.clb,trans=opt.trans,proj=opt.proj,blat=opt.blat,blon=opt.blon,back=opt.back,vmin=opt.vmin,vmax=opt.vmax,ndiv=opt.ndiv,save=opt.save)

def sparse(oplot=None,file=None,ftype=None,xvar=None,yvar=None,zvar=None,cvar=None,xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None,cmin=None,cmax=None,cond=None,merge=False,clb="def",trans=1.,proj=None,blat=None,blon=None,back=None,vmin=None,vmax=None,ndiv=10,save="gui"):

   import sys
   from os import system,path
   from scipy.io.idl import readsav
   import numpy as np
   from myplot import separatenames, definesubplot, errormess, defcolorb, fmtvar, polarinterv, simplinterv, define_proj, readdata, makeplotres
   from mymath import min, max, writeascii
   import matplotlib as mpl
   from matplotlib.pyplot import subplot, figure, plot, scatter, colorbar, show,  title, close, legend, xlabel, ylabel, axis, hist
   from matplotlib.cm import get_cmap
   from mpl_toolkits.mplot3d import Axes3D


   #############################
   ### Load and check data

   #############################
   ### Load and check data

   if file is None:
      print "You must specify at least a file to process with -f."
      exit()
   if xvar is None:
      print "You must specify at least a 1st field with -X."
      exit()
   if proj is not None and yvar is None:
      print "Why did you ask a projection with only one variable?"
      exit()
      
   filename=separatenames(file)
   
   #print 'conditions', cond
   
   print vmax
   print vmin
   
   if cond is None: nslices = 1
   else:                nslices = len(cond)
   numplot   = nslices
   if merge is False: numplot = numplot*len(filename)
   
   print ' '
   if merge is False:
      print nslices, 'condition(s) and', len(filename), 'file(s) without merging ---->', numplot, 'plot(s)'
   else:
      print nslices, 'condition(s) and', len(filename), 'file(s) with merging ---->', numplot, 'plot(s)'
   print ' '

   all_type = [[]]*len(filename)
   all_x  = [[]]*len(filename)
   all_y  = [[]]*len(filename)
   all_z  = [[]]*len(filename)
   all_c  = [[]]*len(filename)
   all_data  = [[]]*len(filename)
   #index  = [[]]*len(filename)


##############################
##################   READ DATA

   for k in range(len(filename)):
   
##### Find file type
     if ftype is None:
      if filename[k].find('.sav') is not -1:
         all_type[k] ='sav'
         print '.sav file detected for', filename[k],'!!'
      elif filename[k].find('.txt') is not -1:
         all_type[k] ='txt'
         #print '.txt file detected for', filename[k],'!!'
      else:
         all_type[k] ='txt'
         print 'no file type detected for', filename[k],'!!, default type is', all_type[k]
     
      if all_type[k] != all_type[0]:
         print 'Two different types were used: ', all_type[k], all_type[0]
         errormess('Not suported')
     else:
        all_type[k] = ftype        

##### Read file       
     if all_type[k] == 'sav':
       print 'reading .sav file', filename[k], '...'
       data = {}
       data = readsav(filename[k], idict=data, python_dict=False) #, verbose=True)
       all_data[k] = data
     elif all_type[k] == 'txt':
       print 'reading .text file', filename[k], '...'
       all_data[k] = np.loadtxt(filename[k])
       
     all_x[k] = readdata(all_data,all_type[k],k,xvar)
     if yvar is not None: all_y[k] = readdata(all_data,all_type[k],k,yvar)
     else:                    all_y[k] = None
     if zvar is not None: all_z[k] = readdata(all_data,all_type[k],k,zvar)
     else:                    all_z[k] = None
     if cvar is not None: all_c[k] = readdata(all_data,all_type[k],k,cvar)
     else:                    all_c[k] = None
       
     if merge is True and k >=1 :
          all_x[0] = np.concatenate((all_x[0],all_x[k]))
          if yvar is not None: all_y[0] = np.concatenate((all_y[0],all_y[k]))
          if zvar is not None: all_z[0] = np.concatenate((all_z[0],all_z[k]))
          if cvar is not None: all_c[0] = np.concatenate((all_c[0],all_c[k]))
       
   
   
##############################
##################   PLOT DATA

## Open a figure and set subplots
   if oplot is not None: fig=oplot
   else: fig = figure()
   subv,subh = definesubplot(numplot,fig)
   palette = get_cmap(name=clb)
   
   
   for nplot in range(numplot):
       
       print ' '
       
###### Find which data and which file for plot nplot
       if merge is False:
          index_s = ((nplot)//len(filename))%nslices
          index_f = ((nplot-1)//nslices)%len(filename)
       else:
          index_s = ((nplot)//1)%nslices
          index_f = ((nplot-1)//nslices)%1
       #print 'nplot,numplot,index_f,index_s', nplot,numplot,index_f,index_s
       
       
###### Select point under conditions defined in -w option
       index = np.isfinite(all_x[index_f])
       if cond is not None:
           zecond = separatenames(cond[index_s])
           #print 'hello', cond[index_s], zecond
           for i in range(len(zecond)):
              zecondi = zecond[i] 
              if zecondi.find('>') != -1:
                 variable = zecondi[0:zecondi.find('>')]
                 value = float(zecondi[zecondi.find('>')+1:len(zecondi)])
                 if merge is True: # ultra moche de reconcatener a chaque fois, mais bon on s'en fout c'est du python
                    zedata = []
                    for k in range(len(filename)):
                       zedata = np.concatenate((zedata,readdata(all_data,all_type[k],k,variable)))
                 else : zedata = readdata(all_data,all_type[k],index_f,variable)
                 #print index.shape,zedata.shape, value
                 index = index*(zedata>value)
                 print 'All points such as', variable, '>', value, 'have been selected'
              elif zecondi.find('<') != -1:
                 variable = zecondi[0:zecondi.find('<')]
                 value = float(zecondi[zecondi.find('<')+1:len(zecondi)])
                 if merge is True:
                    zedata = []
                    for k in range(len(filename)):
                       zedata = np.concatenate((zedata,readdata(all_data,all_type[k],k,variable)))
                 else : zedata = readdata(all_data,all_type[k],index_f,variable)
                 #print index.shape,zedata.shape, value
                 index = index*(zedata<value)
                 print 'All points such as', variable, '<', value, 'have been selected'
              else:
                 print ''
                 print 'I do not understand that condition :', zecondi
                 errormess('')
       
       
       if np.sum(index) == 0:
           print '***********  WARNING  ***********'
           print '***********  NO DATA  ***********'
           errormess('')
       else:
           print np.sum(index),'points have been selected among', len(all_x[index_f]), \
           'that is to say %2.2f' %(float(100*np.sum(index))/float(len(all_x[index_f]))), '%.'
       
          
       #print 'numplot', numplot
       changesubplot = (numplot > 1) #and (len(what_I_plot.shape) != 1)
       if changesubplot: subplot(subv,subh,nplot+1)
       #print 'subv,subh,nplot', subv,subh,nplot
       
###### Projection
       if proj is not None: # Nota : NEVER TRY TO DO A MESHGRID ON SPARSE DATA. WAY TOO MUCH POINTS.
          wlon = [min(all_x[index_f][index]),max(all_x[index_f][index])]
          wlat = [min(all_y[index_f][index]),max(all_y[index_f][index])]
          m = define_proj(proj,wlon,wlat,back=back,blat=blat,blon=blon)
          x, y = m(all_x[index_f][index],all_y[index_f][index])
       else:
          x = all_x[index_f][index]
          if yvar is not None: y = all_y[index_f][index]                              
       
###### Plot: 1D histogram
       if yvar is None:
          hist(x,bins=ndiv,histtype='step',linewidth=2)
          if save == 'txt':
              print 'saving file profile'+str(nplot+1)+'.txt'
              writeascii(np.transpose(x),'profile'+str(nplot+1)+'.txt')
###### Plot: 2D cloud
       elif zvar is None and cvar is None :
          plot(x,y,'.b')
          if save == 'txt':
               print 'saving file profile'+str(nplot+1)+'.txt'
               writeascii(np.transpose(np.array([x,y])),'profile'+str(nplot+1)+'.txt')
###### Plot: 2D cloud with color
       elif zvar is None and cvar is not None :
          if save == 'txt':
               print 'saving file profile'+str(nplot+1)+'.txt'
               writeascii(np.transpose(np.array([x,y,all_c[index_f][index]])),'profile'+str(nplot+1)+'.txt')
          scatter(x,y,c=all_c[index_f][index],\
          marker='o', edgecolor='None',cmap = palette, alpha=trans, vmin = vmin,vmax = vmax)
###### Plot: 3D cloud
       elif zvar is not None and cvar is None :
          ax = fig.add_subplot(subv,subh,nplot+1, projection='3d')
          ax.plot(x,y,all_z[index_f][index],'.b')
###### Plot: 3D cloud with color
       else :
          ax = fig.add_subplot(subv,subh,nplot+1, projection='3d')
          ax.scatter(x,y,all_z[index_f][index],c=all_c[index_f][index],\
          marker='o', edgecolor='None', cmap = palette, alpha=trans,vmin = vmin,vmax = vmax)
   
###### Colorbar: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps?action=AttachFile&do=get&target=colormaps3.png
       if clb != 'nobar' and all_c[index_f] is not None and all_z[index_f] is None: # pourquoi la colorbar marche pas en 3d?
          colorbar( fraction=0.05,pad=0.03,format='%0.2f',\
          extend='neither',spacing='proportional' )
        
###### Plot limits (options --xmin, --xmax, etc ..) 
       if proj is None:
          xlabel(xvar)
          ylabel(yvar)
          if zvar is None :
             if xmin is not None: mpl.pyplot.xlim(xmin=xmin)
             else:                    mpl.pyplot.xlim(xmin=min(all_x[index_f][index]))
             if xmax is not None: mpl.pyplot.xlim(xmax=xmax)
             else:                    mpl.pyplot.xlim(xmax=max(all_x[index_f][index]))
             if yvar is not None:
                if ymin is not None: mpl.pyplot.ylim(ymin=ymin)
                else:                    mpl.pyplot.ylim(ymin=min(all_y[index_f][index]))
                if ymax is not None: mpl.pyplot.ylim(ymax=ymax)
                else:                    mpl.pyplot.ylim(ymax=max(all_y[index_f][index]))

          
       if cond is not None:
          title(cond[index_s])
       else:
          title('all point selected')
       
   zeplot = "output"
   if save in ['png','eps','svg','pdf']:     makeplotres(zeplot,ext=save)#,res=resolution,pad_inches_value=pad_inches_value,disp=False,ext=save)
   elif save == 'gui':                       show()
   elif save == 'return':                    return mpl.pyplot
   else:                                         print "INFO: save mode not supported. using gui instead." ; show()
   print ''

   command = ""
   for arg in sys.argv: command = command + arg + ' '
   name = zeplot
   f = open(name+'.sh', 'w')
   f.write(command)

   


   
