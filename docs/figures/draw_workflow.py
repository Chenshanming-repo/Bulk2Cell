"""Render the bulk2cell README workflow as editable SVG and high-resolution PNG."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

OUT = Path(__file__).resolve().parent
plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(16, 10), facecolor='white')
ax.set(xlim=(0, 16), ylim=(0, 10))
ax.axis('off')
black, blue, orange, green, purple = '#303842', '#4077AD', '#D99020', '#19945A', '#8C4799'

def text(x, y, s, size=12, weight='normal', color=black, ha='center'):
    ax.text(x, y, s, fontsize=size, weight=weight, color=color,
            ha=ha, va='center', linespacing=1.45)

def route(points, color, arrow=True, lw=5):
    xs, ys = zip(*points)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle='round', zorder=1)
    if arrow:
        a, b = points[-2:]
        dx, dy = b[0]-a[0], b[1]-a[1]
        start = (b[0]-dx*.15, b[1]-dy*.15)
        ax.add_patch(FancyArrowPatch(start, b, arrowstyle='-|>', mutation_scale=19,
                                    lw=0, color=color, zorder=2))

def node(x, y, label, below=True):
    ax.scatter([x], [y], s=180, c=purple, edgecolors='white', linewidths=2, zorder=3)
    text(x, y+(-.49 if below else .49), label, 12, 'bold')

def input_label(x, y, title, detail):
    text(x, y, title, 13, 'bold', '#B52D50')
    text(x, y-.32, detail, 10, color='#65717E')

text(.5, 9.55, 'bulk2cell', 29, 'bold', ha='left')
text(.5, 9.08, 'Bulk long-read-assisted isoform quantification at single-cell resolution', 14, ha='left')

input_label(1.7, 8.3, 'PacBio bulk reads', 'HiFi + primers / FLNC BAMs')
route([(3.1, 8.25), (14, 8.25)], black)
node(4, 8.25, 'lima + refine\n(HiFi input)')
node(6.5, 8.25, 'Iso-Seq cluster\n+ genome alignment')
node(9, 8.25, 'Collapse\n+ Pigeon filter')
node(12, 8.25, 'Union transcript catalog\nInput GTF + Iso-Seq')
route([(14, 8.25), (14.8, 8.25), (14.8, 4.35), (12.5, 4.35)], black)
text(15.12, 6.3, 'GTF\n+\nFASTA', 10, color=black)

input_label(1.7, 6.35, 'Original reference', 'Genome FASTA + annotation GTF')
route([(3.1, 6.3), (4.3, 6.3), (4.3, 5.65), (12.5, 5.65)], blue)
node(4.3, 6.3, 'Cell Ranger mkref', below=False)
node(7.2, 5.65, 'Cell Ranger count')
node(10, 5.65, 'Tagged BAM\n+ called barcodes', below=False)
node(12.5, 5.65, 'Called-cell\nRNA pseudobulk', below=False)
input_label(7.2, 6.8, 'Illumina 10x reads', 'Paired FASTQs + chemistry')
route([(7.2, 6.35), (7.2, 5.75)], blue, lw=3)
text(1.7, 5.45, 'Reference also supports\nIso-Seq alignment, Pigeon\nand union construction.', 10, color='#65717E')

route([(12.5, 5.4), (12.5, 4.35)], orange)
node(12.5, 4.35, 'Salmon index + quant', below=True)
route([(12.3, 4.35), (7.2, 4.35), (7.2, 2.8)], orange)
text(8.35, 4.64, 'Transcript abundance prior', 11, color=orange)
route([(10, 5.42), (10, 3.3), (9.2, 3.3), (9.2, 2.8)], blue)
text(10.65, 3.65, 'Cell / UMI\nevidence', 10, color=blue)
route([(14.8, 4.35), (14.8, 3.15), (11.2, 3.15), (11.2, 2.8)], black)
text(13.15, 3.38, 'Reference + Iso-Seq models', 10)

route([(11.2, 2.8), (2.7, 2.8)], green)
node(11.2, 2.8, 'Molecule–isoform\ncompatibility')
node(7.2, 2.8, 'Capture model\n+ group EM')
node(3.2, 2.8, 'Cell × isoform / group\ncount matrices + QC')
text(7.3, 1.65, 'bulk2cell quantify-full', 15, 'bold', green)

# Small matrix motif represents outputs without inventing biological results.
for row in range(3):
    for col in range(5):
        ax.add_patch(Rectangle((.65+col*.23, 2.75+row*.23), .2, .2,
            facecolor=['#E0EFE7','#94C8AB','#47996D'][(row+col)%3], edgecolor='none'))

for x, c, label in [(1.0, black, 'Transcript catalog'), (4.9, blue, 'Illumina processing'),
                     (9.0, orange, 'Abundance prior'), (12.4, green, 'Quantification')]:
    ax.plot([x, x+.45], [.85, .85], lw=5, color=c)
    text(x+.6, .85, label, 10, ha='left')
text(8, .3, 'Cell Ranger uses only the original reference and Illumina reads; it has no Iso-Seq dependency.', 11)
fig.subplots_adjust(left=.01, right=.99, top=.99, bottom=.01)
for suffix in ('svg', 'png'):
    fig.savefig(OUT/f'bulk2cell_workflow.{suffix}', dpi=180, facecolor='white')
plt.close(fig)
