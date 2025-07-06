#=============================================================
#
# Import Seaborn to make a pretty plot.
#
# Seaborn is excellent for controlling decorations, line widths
# color palettes, etc.
#
#   - Define background for the plot
#   - Define if a grid is used (background)
#   - Define font families
#   - Define line widths
#
# Ideally, this should be put in a separate file and imported, 
# but is included here as an example.
#=============================================================

import seaborn as sns
from matplotlib import rc


sns.set()
sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.0})

#---- This allows the use of LaTeX + the use sans-serif fonts also for tick labels:

rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')

# There are 5 presents for background: darkgrid, whitegrid, dark, white, and ticks
# Define how ticks are placed and define font families

sns.set_style("ticks")
sns.set_style("whitegrid", 
 {'axes.edgecolor': 'black',
 'axes.grid': True,
 'axes.axisbelow': True,
 'axes.labelcolor': '.15',
 'grid.color': '0.9',
 'grid.linestyle': '-',
 'xtick.direction': 'in', 
 'ytick.direction': 'in',
 'xtick.bottom': True,
 'xtick.top': True,
 'ytick.left': True,
 'ytick.right': True, 
 'font.family': ['sans-serif'],
 'font.sans-serif': [
  'Liberation Sans',
  'Bitstream Vera Sans',
  'sans-serif'],})
