import sys,os,commands
import  os
import re
import sys

# environment
USE_CULA=False
USE_CIMG=False
USE_OPENCV=False

#print "----------------------------------------------"
Decider('MD5-timestamp')
LOCAL_ROOT="/opt/local"
DHCAL_ROOT=os.path.abspath("../..")

# includes
INCLUDES=['.','./include',"/usr/include","/usr/local/include","/usr/X11R6/include","/usr/include/jsoncpp"]

#MYSQL and sqlite
INCLUDES.append("/opt/zdaq/include")
INCLUDES.append("/opt/lydaq/drivers/common/include")
INCLUDES.append("/usr/include/libxml2")
#python
INCLUDES.append(commands.getoutput("python -c 'import distutils.sysconfig as conf; print conf.get_python_inc()'"))

#flags
CPPFLAGS=commands.getoutput("$ROOTSYS/bin/root-config --cflags").split(" ") + ["-g","-DLINUX", "-DREENTRANT" ,"-Dlinux", "-DLITTLE_ENDIAN__", "-Dx86","-DDAQ_VERSION_2","-Wno-deprecated","-std=c++11"] 
if (USE_OPENCV):
	CPPFLAGS.append(commands.getoutput("pkg-config --cflags opencv").split(" "));
	CPPFLAGS.append("-DUSE_OPENCV")
if (USE_CIMG):
	CPPFLAGS.append("-Dcimg_use_vt100")
	CPPFLAGS.append("-DUSE_CIMG")
if (USE_CULA):
	CPPFLAGS.append("-DUSE_CULA")
#Library ROOT + some of XDAQ + DB 
ROOT_LIBS=[lib[2:] for lib in filter(lambda x: (x[:2]=="-l"), commands.getoutput("$ROOTSYS/bin/root-config --libs --ldflags ").split(" "))
]
ROOT_LIBS.append('XMLIO')
ROOT_LIBPATH=[lib[2:] for lib in filter(lambda x: (x[:2]=="-L"), commands.getoutput("$ROOTSYS/bin/root-config --libs  ").split(" "))]
# Summary

LIBRARIES=ROOT_LIBS+['pthread','jsoncpp','zdaq','log4cxx']



#Library path XDAQ,DHCAL and ROOT + Python
LIBRARY_PATHS=ROOT_LIBPATH+["/usr/lib","/usr/local/lib","/opt/dhcal/lib","./lib"]
LIBRARY_PATHS.append(commands.getoutput("python -c 'import distutils.sysconfig as conf; print conf.PREFIX'")+"/lib")

#link flags

LDFLAGS=commands.getoutput("$ROOTSYS/bin/root-config --ldflags").split(" ")
LDFLAGS.append(['-fPIC','-dynamiclib','-rdynamic'])
		      
# Create the Environment
env = Environment(CPPPATH=INCLUDES,CPPFLAGS=CPPFLAGS,LINKFLAGS=LDFLAGS, LIBS=LIBRARIES,LIBPATH=LIBRARY_PATHS)

#print "CC is:",env.subst('$CPPPATH')




#Shared library
lmrb=env.SharedLibrary("#lib/lmrb",source=['src/jsonGeo.C','src/DCHistogramHandler.C','src/TdcAnalyzer.C','src/tdcrb.cxx','src/RecoAnalyzer.C','src/RbServer.cxx'])
LIBRARIESPLUG=["lmrb","boost_system"]+LIBRARIES
LIBRARY_PATHS.append(".")

env1 = Environment(CPPPATH=INCLUDES,CPPFLAGS=CPPFLAGS,LINKFLAGS=LDFLAGS, LIBS=LIBRARIESPLUG,LIBPATH=LIBRARY_PATHS)






EXE_LIBPATH=LIBRARY_PATHS

EXE_LIBS=LIBRARIESPLUG


tdcr=env1.Program("bin/tdcr",source="src/tdcr.cc",LIBPATH=EXE_LIBPATH,LIBS=EXE_LIBS)

rbd=env1.Program("bin/rb_daemon",source="src/RbDaemon.cc",LIBPATH=EXE_LIBPATH,LIBS=EXE_LIBS)





Default([lmrb,tdcr,rbd])




#env1.Install("../../opt/dhcal/lib",[lmreco,rawAnalysis,trackAnalysis,rootProcessor,combAnalysis])
#env1.Alias('install', ["../../opt/dhcal/lib"])
