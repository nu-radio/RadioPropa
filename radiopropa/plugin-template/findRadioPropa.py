# Returns the path to the radiopropa swig_interface, install_prefix, or -1 if
# radiopropa is not found
import sys
try:
  import radiopropa
except ImportError:
  sys.exit(-1)

if sys.argv[1] == 'swig_interface':
  sys.stdout.write(radiopropa.getDataPath('swig_interface'))
elif sys.argv[1] == 'install_prefix':
  sys.stdout.write(radiopropa.getInstallPrefix())
