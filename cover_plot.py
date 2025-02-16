from statistics import mean
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import argparse


parser = argparse.ArgumentParser(
    description="If you used samtools mpileup, you can plot it")
parser.add_argument("--strain", type=str, required=False, default='',
                    help="Name of strain that wil be printed in title")
parser.add_argument("--chr", type=str, required=False,
                    default='chromosome',
                    help="Name of chromosome")
parser.add_argument("--cntrchr", type=str, required=False, default='1',
                    help="Name of control chromosome, on which cover "
                    "will be normolized. If it is unused, plot will be "
                    "whith real cover.")
parser.add_argument("--chrcov", type=str, required=True,
                    help="Path to chromosom cover")
parser.add_argument("--cntrcov", type=str, required=False, default='',
                    help="Path to control chromosome cover")
parser.add_argument("--window", type=int, required=False, default=1000,
                    help="Size of window for moving average, "
                    "default=1000")
parser.add_argument("--title", type=str, required=False, default=None,
                    help="Customize title")
parser.add_argument("--xlabel", type=str, required=False, default=None,
                    help="Customize xlabel")
parser.add_argument("--legend", type=str, required=False, default=[], nargs="+",
                    help="Customize legend")
parser.add_argument("--lang", type=str, required=True,
                    help="RU or EN")
parser.add_argument("--out", type=str, required=False, default="./",
                    help='Path to output default="./"')
parser.add_argument("--i", type=str, required=False, default="zero",
                    help="Any ID for plot")

args = parser.parse_args()

WINDOW = args.window

def moving_average(Y, window):
    averaged = []
    L = len(Y)
    
    sma_prev = sum(Y[:window])/window
    averaged.append(sma_prev)
    
    for i in range(window, L):
        sma_i = sma_prev + (Y[i] - Y[i-window])/window
        averaged.append(sma_i)
        sma_prev = sma_i
    
    return averaged

#average cover on cntr chr
if args.cntrcov == '':
    average_cover = 1
else:
    covers = []
    with open(args.cntrcov) as inp:
        for line in inp:
                line = line.split("\t")
                covers.append(int(line[3]))
        average_cover = mean(covers)

#cover on chr
with open(args.chrcov) as inp:
    poses =[]
    covers=[]
    for line in inp:
            line = line.split("\t")
            poses.append(int(line[1]))
            covers.append(int(line[3])/average_cover)
    

    
#window
new_poses = poses[(WINDOW-1):]
new_covers = moving_average(covers, WINDOW)

plt.plot(poses, covers, color="#BF9315")    
plt.plot(new_poses, new_covers, color="#7B018C")   
plt.plot(range(0, max(poses)), [1 for x in range(0, max(poses))],
         color="#0F8C69")


if args.lang == "EN":    
    if args.cntrchr:
        plt.title(f"Coverage for {args.strain} on {args.chr}\n"
                  f"normalized to {args.cntrchr}\n", fontweight="bold")
    else:
        plt.title(f"Coverage for {args.strain} on {args.chr}\n",
                  fontweight="bold")
    plt.xlabel(f"Position in {args.chr}")
    plt.legend(['All positions',
                f'In {WINDOW} window',
                f'Average cover on\n'
                f'{args.cntrchr}: {average_cover:.02f}'
                ], loc='center left', bbox_to_anchor=(1, 0.5),
               )


if args.lang == "RU":
    if args.cntrchr:
        plt.title(f"Покрытие {args.chr} для {args.strain} \n"
                  f"нормированное на {args.cntrchr}\n", fontweight="bold")
    else:
        plt.title(f"Покрытие {args.chr} для {args.strain}\n",
                  fontweight="bold")
    plt.xlabel(f"Позиция в {args.chr}")
    plt.legend(['Все позиции',
                f'В окне {WINDOW}',
                f'Среднее покрытие\n'
                f'{args.cntrchr}: {average_cover:.02f}'
                ], loc='center left', bbox_to_anchor=(1, 0.5),
               )


if args.title:
    plt.title(f"{args.title}", fontweight="bold")
if args.xlabel:
    plt.xlabel(f"{args.xlabel}")
if args.legend:
    plt.legend(args.legend, loc='center left', bbox_to_anchor=(1, 0.5))       

plt.savefig(f"{args.out}/{args.strain}_plotty_{args.lang}_{args.i}.svg", bbox_inches="tight")
plt.savefig(f"{args.out}/{args.strain}_plotty_{args.lang}_{args.i}.png", bbox_inches="tight")
