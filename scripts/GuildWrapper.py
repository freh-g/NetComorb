import subprocess
import argparse
import os

parser=argparse.ArgumentParser(description='parse the HIPPIE PPI dataset and DisGeNet database')
parser.add_argument('-a','--algorithm',type=str,help='which algorithm to use')
parser.add_argument('-d','--disease',type=str,help='disease Name useful for multiple phenotypes')
args=parser.parse_args()



def run_NetScore(PathGuild,PathToNode,PathToInteractome,PathOut,r=3,i=2):
    
    OutDir='../outputs/NetScore/'
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    
    subprocess.run('bash -c "source activate guildmaster && %s./guild_x64 -s s -n %s -e %s -o %s -r %s -i %s"' % (PathGuild,PathToNode,PathToInteractome,PathOut,r,i),shell=True)




def create_randomized_networks(ScriptPath,PathInteractome,n=100):
    """This script creates randomized interactomes in the same folder of the interactome and then move
    them in a folder called randomized interactomes choosen by the user"""
    
    
    ##build the randomized interactome networks
    subprocess.call('bash -c "source activate guildmaster && python %s %s %s"' %(ScriptPath,PathInteractome,n),
                   shell=True)
    subprocess.run('mv %s.* %s'%(PathInteractome,PathOut),shell=True)


def run_NetZcore(PathGuild,PathNode,PathInteractome,PathOut,PathInteractomes,i=5):
    OutDir='../outputs/NetZcore/'
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    subprocess.run('bash -c "source activate guildmaster && %s./guild_x64 -s z -n %s -e %s -o %s -i %s -d %s -x 100"' % (PathGuild,PathNode,PathInteractome,PathOut,i,PathInteractomes),shell=True)


if __name__=='__main__':
    
    if args.algorithm=='s':
        
        run_NetScore(
            '../guild/',
            '../inputs/%s_nodes.sif'%(args.disease),
            '../inputs/hippie_interactome.sif',
            '../outputs/NetScore/%s_output_NetScore'%(args.disease)
            )
    
    elif args.algorithm=='z':
        
        PathOut='../inputs/RandomizedNetworks/'
        if not os.path.exists(PathOut):
            os.makedirs(PathOut)
        if len(os.listdir(PathOut))==0:
            create_randomized_networks(

                '../guild/src/create_random_networks_for_netzcore.py',
                '../inputs/hippie_interactome.sif',
                '100'

                )
        
        run_NetZcore(
            
            '../guild/',
            '../inputs/%s_nodes.sif'%(args.disease),
            '../inputs/hippie_interactome.sif',
            '../outputs/NetZcore/%s_output_NetZcore'%(args.disease),
            '../inputs/randomizedNetworks/hippie_interactome.sif.'

        )
        
