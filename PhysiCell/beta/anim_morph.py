#
# anim_morph.py:  render/animate PhysiCell .mat (cell) files, using left/right arrows on keyboard
#
# Keyboard arrows: right/left arrows will single step forward/backward; up/down will increment/decrement step size
#
# Dependencies include matplotlib and numpy. We recommend installing the Anaconda Python3 distribution.
#
__author__ = "Randy Heiland"

import sys
import os
import math
from pyMCDS_optional_meshes import pyMCDS
join_our_list = "(Join/ask questions at https://groups.google.com/forum/#!forum/physicell-users)\n"
try:
  import matplotlib
  import matplotlib.colors as mplc
  from matplotlib.patches import Circle, Ellipse, Rectangle
  from matplotlib.collections import PatchCollection
except:
  print("\n---Error: cannot import matplotlib")
  print("---Try: python -m pip install matplotlib")
  print(join_our_list)
#  print("---Consider installing Anaconda's Python 3 distribution.\n")
  raise
try:
  import numpy as np  # if mpl was installed, numpy should have been too.
except:
  print("\n---Error: cannot import numpy")
  print("---Try: python -m pip install numpy\n")
  print(join_our_list)
  raise
from collections import deque
try:
  # apparently we need mpl's Qt backend to do keypresses 
#  matplotlib.use("Qt5Agg")
  matplotlib.use("TkAgg")
  import matplotlib.pyplot as plt
except:
  print("\n---Error: cannot use matplotlib's TkAgg backend")
  print(join_our_list)
#  print("Consider installing Anaconda's Python 3 distribution.")
  raise


current_idx = 0
print("# args=",len(sys.argv)-1)

#for idx in range(len(sys.argv)):
use_defaults = True
show_nucleus = 0
current_idx = 0
# axes_min = 0.0
# axes_max = 1000  # but overridden by "width" attribute in .svg
axes_min = -200.0
axes_max = 200  # but overridden by "width" attribute in .svg
if (len(sys.argv) == 5):
  use_defaults = False
  kdx = 1
  show_nucleus = int(sys.argv[kdx])
  kdx += 1
  current_idx = int(sys.argv[kdx])
  kdx += 1
  axes_min = float(sys.argv[kdx])
  kdx += 1
  axes_max = float(sys.argv[kdx])
elif (len(sys.argv) != 1):
  print("Please provide either no args or 4 args:")
  usage_str = "show_nucleus start_index axes_min axes_max"
  print(usage_str)
  print("e.g.,")
  eg_str = "%s 0 0 0 2000" % (sys.argv[0])
  print(eg_str)
  sys.exit(1)

#"""
print("show_nucleus=",show_nucleus)
print("current_idx=",current_idx)
print("axes_min=",axes_min)
print("axes_max=",axes_max)
#"""

"""
if (len(sys.argv) > 1):
   current_idx = int(sys.argv[1])
if (len(sys.argv) > 2):
   axes_min = float(sys.argv[2])
   axes_max = float(sys.argv[3])

if (len(sys.argv) > 4):
  usage_str = "[<start_index> [<axes_min axes_max>]]"
  print(usage_str)
  print("e.g.,")
  eg_str = "%s 1 10 700 1300" % (sys.argv[0])
  print(eg_str)
  sys.exit(1)
"""

print("current_idx=",current_idx)

#d={}   # dictionary to hold all (x,y) positions of cells

""" 
--- for example ---
In [141]: d['cell1599'][0:3]
Out[141]: 
array([[ 4900.  ,  4900.  ],
       [ 4934.17,  4487.91],
       [ 4960.75,  4148.02]])
"""

fig = plt.figure(figsize=(7,7))
ax = fig.gca()
#ax.set_aspect("equal")


#plt.ion()

time_delay = 0.1

count = -1
#while True:

#-----------------------------------------------------
#def ellipses(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
def ellipses(x, y, width, height, angle, c='b', vmin=None, vmax=None, **kwargs):
    """
    See https://gist.github.com/syrte/592a062c562cd2a98a83 

    Make a scatter plot of circles. 
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    # zipped = np.broadcast(x, y, width, height, angle, s)
    zipped = np.broadcast(x, y, width, height, angle)
    # patches = [Ellipse((x_, y_), width_, height_, angle_, s_)
    patches = [Ellipse((x_, y_), width_, height_, angle_)
               for x_, y_, width_, height_, angle_  in zipped]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    plt.draw_if_interactive()
    if c is not None:
        plt.sci(collection)
    return collection

#-----------------------------------------------------
def plot_cells():
  global current_idx, axes_max

  fname = "output%08d.xml" % current_idx
  if (os.path.isfile(fname) == False):
    print("File does not exist: ",fname)
    return
  mcds = pyMCDS(fname,'.')  
  tmins = mcds.get_time()
  print('time (mins)=',tmins)
  title_str = "Current time: " + str(tmins)

  num_cells = len(mcds.data['discrete_cells']['ID'])
  print('num_cells = ',num_cells)

  xvals = mcds.data['discrete_cells']['position_x']
  print("xvals= ",mcds.data['discrete_cells']['position_x'])
  yvals = mcds.data['discrete_cells']['position_y']

  axis_a = mcds.data['discrete_cells']['axis_a']
  axis_b = mcds.data['discrete_cells']['axis_b']
  axis_c = mcds.data['discrete_cells']['axis_c']

  xlist = deque()
  ylist = deque()
  wlist = deque()
  hlist = deque()
  alist = deque()
  rgb_list = deque()

  height = 10
  height_del = 5  # vary just to see different sizes
  width = 30
  width_del = 10  # vary just to see different sizes
  angle = tmins * 45


#  print('--- child.tag, child.attrib ---')
  numChildren = 0
  for icell in range(num_cells):
      xval = xvals[icell]
      yval = yvals[icell]

      s = 'red'
      rgb_tuple = mplc.to_rgb(mplc.cnames[s])  # a tuple
      rgb = [x for x in rgb_tuple]

      xlist.append(xval)
      ylist.append(yval)

    #   width += width_del
      a_val = axis_a[icell]
      wlist.append(a_val)    # "width" of each ellipse/cell

      height += height_del
      hlist.append(height)    # "height" of each ellipse/cell
      alist.append(angle + icell*45)
      rgb_list.append(rgb)

  xvals = np.array(xlist)
  yvals = np.array(ylist)
  widths = np.array(wlist)
  heights = np.array(hlist)
  angles = np.array(alist)
  rgbs =  np.array(rgb_list)

  plt.cla()
  title_str += " (" + str(num_cells) + " agents)"
  plt.title(title_str)
  plt.xlim(axes_min,axes_max)
  plt.ylim(axes_min,axes_max)

  ellipses(xvals,yvals, widths, heights, angles, color=rgbs)
  plt.pause(time_delay)

step_value = 1
def press(event):
  global current_idx, step_value
#    print('press', event.key)
  sys.stdout.flush()
  if event.key == 'escape':
    sys.exit(1)
  elif event.key == 'h':  # help
    print('esc: quit')
    print('right arrow: increment by step_value')
    print('left arrow:  decrement by step_value')
    print('up arrow:   increment step_value by 1')
    print('down arrow: decrement step_value by 1')
    print('0: reset to 0th frame')
    print('h: help')
  elif event.key == 'left':  # left arrow key
#    print('go backwards')
#    fig.canvas.draw()
    current_idx -= step_value
    if (current_idx < 0):
      current_idx = 0
    plot_cells()
  elif event.key == 'right':  # right arrow key
#        print('go forwards')
#        fig.canvas.draw()
    current_idx += step_value
    plot_cells()
  elif event.key == 'up':  # up arrow key
    step_value += 1
    print('step_value=',step_value)
  elif event.key == 'down':  # down arrow key
    step_value -= 1
    if (step_value <= 0):
      step_value = 1
    print('step_value=',step_value)
  elif event.key == '0':  # reset to 0th frame/file
    current_idx = 0
    plot_cells()
  else:
    print('press', event.key)


#for current_idx in range(40):
#  fname = "snapshot%08d.svg" % current_idx
#  plot_cells(fname)
plot_cells()
print("\nNOTE: click in plot window to give it focus before using keys.")

fig.canvas.mpl_connect('key_press_event', press)

# keep last plot displayed
#plt.ioff()
plt.show()
