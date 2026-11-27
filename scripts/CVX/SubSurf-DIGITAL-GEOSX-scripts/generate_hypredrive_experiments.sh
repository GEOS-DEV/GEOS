#!/bin/bash

# $0 "N_cores ..." "N_nodes"

source ${HOME}/cs691/utils/bin/utils_bash.sh

if [[ $# < 3 || -z $INITRC ]]; then
    RCerror $0 10 "Usage : INITRC='/data/saet/software/.init/hypredrive/XXXX.rc ' [WHERE='Slurm'] [N_y=][N_z=] $0 \"N_cores ...\" \"N_nodes ...\" \"N_elems_core ...\" "
fi

: ${V:="0"}; 

[[ $V -ge 10 ]] && set -x
: ${CLIMPIOPTS:=" -mca coll_hcoll_enable 0 --report-bindings --bind-to core --oversubscribe "};
export CLIMPIOPTS

: ${USEGPUS:="0"};
export USEGPUS

: ${stencil:=7}
export stencil

: ${n_x:=1}
: ${N_y:=100}
: ${N_z:=100}
export n_x N_y N_z

: ${WHERE:="Slurm"};

if [[ ${WHERE} == "Slurm" ]]; then
    sbatch="1";
else
    sbatch="0";
fi
export WHERE

export N_cores="$1"
export N_nodes="$2"
export N_elems_core="$3"
 
if [[ ${WHERE} == "localhost" ]]; then
    N_nodes=1 ;
fi
    
export thr_rank=" 1 "

typeset -A SKU
SKU["1"]="a100"
SKU["2"]="a100"
SKU["3"]="a100"
SKU["4"]="a100"
SKU["5"]="a100"
SKU["6"]="a100"
SKU["7"]="a100"
SKU["8"]="a100"
SKU["120"]="hbv3"
SKU["176"]="hbv4"
SKU["352"]="hbv5"

export SKU

typeset -A partition
partition["1"]="geos";
partition["2"]="geos";
partition["3"]="geos";
partition["4"]="geos";
partition["5"]="geos";
partition["6"]="geos";
partition["7"]="geos";
partition["8"]="geos";
partition["120"]="geos";
partition["176"]="geos";
partition["352"]="hbv5test";

export partition

echo "#! /bin/bash "
echoc0 "=======================================================================================#"
echoc0 "N_cores=\"$N_cores\" "
echoc0 "N_nodes=\"$N_nodes\" "
echoc0 "N_elems_core=\"$(echo $N_elems_core)\" "
echoc0 "=======================================================================================#"

echo "export INITRC=$INITRC"
echo "export CLIMPIOPTS=\"$CLIMPIOPTS\""
echo "export PINDOM=\"numa\""

for n_cores in $N_cores; do
    export n_2cores=$(echo "$n_cores*1/2" | bc )
    export n_3cores=$(echo "$n_cores*3/4" | bc )
    [[ $V -ge "1" ]] && echoc0 "n_cores=$n_cores  n_2cores=$n_2cores  n_3cores=$n_3cores"
    for n_elems_core in $N_elems_core; do 
	[[ $V -ge "1" ]] && echoc0 "n_elems_core=$n_elems_core"
	for n_nodes in $N_nodes; do
	    for n_y in $N_y; do 
		for n_z in $N_z; do 
		    export n_yz=$(echo "$n_y*$n_z"| bc )
		    N_elems_node=$(echo "$n_cores*$n_elems_core" | bc ) ;
		    N_total_cores=$(echo "$n_cores*$n_nodes" | bc ) ;
		    N_total_elems=$(echo "$N_total_cores*$n_elems_core" | bc ) ;
		    N_elems_2node=$(echo "$n_2cores*$n_elems_core" | bc ) ;
		    N_total_2cores=$(echo "$n_2cores*$n_nodes" | bc ) ;
		    N_total_2elems=$(echo "$N_total_2cores*$n_elems_core" | bc ) ;
		    N_elems_3node=$(echo "$n_3cores*$n_elems_core" | bc ) ;
		    N_total_3cores=$(echo "$n_3cores*$n_nodes" | bc ) ;
		    N_total_3elems=$(echo "$N_total_3cores*$n_elems_core" | bc ) ;
		    N_elems_cube=$(echo "e(l($N_total_elems)/3)" | bc -l | tr '.' ' ' | gawk '{print $1}')

		    [[ $V > "1" ]] && { echoc0 "N_total_cores= $N_total_cores" ;
					echoc0 "N_total_2cores= $N_total_2cores" ;
					echoc0 "N_total_3cores= $N_total_3cores" ;
					echoc0 "N_total_elems= $N_total_elems" ;
					echoc0 "N_elems_cube= $N_elems_cube"
		    }
		    n_x=$(echo "$N_total_elems $n_yz" | gawk '{printf "%d", $1/$2 }'); [[ $nx -lt 1 ]] && n_x=1 ; 
		    [[ $V > "1" ]] && echoc0 " n_x = $n_x " ;
		    export n_act_elems=$(echo "$n_x * $n_y * $n_z "| bc )
		    eid="${N_total_cores}--${n_x}x${n_y}x${n_z}" 

		    P_z=$( echo $n_z $n_elems_core |
			       gawk '{n_z=$1; n_elems_core=$2; P_z = int(n_z / n_elems_core); if (P_z > n_z) print n_z; else print (P_z * n_elems_core < n_z) ? (P_z + 1) : P_z }' )
		    P_y=$( echo $n_y $n_elems_core |
			       gawk '{n_y=$1; n_elems_core=$2; P_y = int(n_y / n_elems_core); if (P_y > n_y) print n_y; else print (P_y * n_elems_core < n_y) ? (P_y + 1) : P_y }' )
		    P_x=$( echo $n_x $N_total_cores $P_y $P_z |
			       gawk -v V=$V '
			       { n_x=$1; N_total_cores=$2; P_y=$3; P_z=$4; P_x = int( N_total_cores / (P_y * P_z) ); 
    			       	 if (P_x > n_x) 
				       P_x = n_x; 
				    else
				       P_x = (P_x * n_elems_core < n_x) ? (P_x + 1) : P_x;
		       	       	 print (P_x > N_total_cores ) ? N_total_cores : P_x ; 
			       }' )

		    N_ranks=$(echo "$P_x * $P_y * $P_z" | bc )

		    appargs=" -P $P_x $P_y $P_z -n ${n_x} ${n_y} ${n_z} -s $stencil # -N $n_nodes ${SKU[${n_cores}]}"
		    echo "# Data Matrix: ${n_act_elems}=$n_x x $n_y x $n_z ; Rank Matrix: $P_x x $P_y x $P_z; N_ranks: $N_ranks ; n_nodes: $n_nodes ; N_total_elems: $N_total_elems ; n_elems_cube: $n_elems_cube $eid"
		    if [[ $P_x -gt $n_x ||  $P_y -gt $n_y ||  $P_z -gt $n_z ]]; then
			echoc0 " Rejected P( $P_x $P_y $P_z ) n( ${n_x} x ${n_y} x ${n_z} )"
		    else
			if [[ ${sbatch} == "1" ]]; then 
			    sbatchcli="sbatch -p ${partition[${n_cores}]} -N $n_nodes --constraint=${SKU[${n_cores}]} --exclusive"
			else
			    sbatchcli=""
			fi
			    
			runcli="USEGPUS=$USEGPUS EXPD=stencil-$stencil EID=$eid ${sbatchcli}  ./run_all-MPI-CUDA.sh $WHERE ' ' laplacian \" $N_ranks \" \" 1 \" $appargs " 
			echo $runcli
		    fi

		done
	    done
	done
    done

done
	
[[ $V -ge 10 ]] && set +x
	
# V=1 EXPD="stencil-7" EID="700-400x200x80" sbatch -p geos --constraint=hbv3 -N 6  --exclusive  ./run_all-MPI-CUDA.sh Slurm " " laplacian " 700 " " 1 "  -P 20 7 5 -n 400 200 80  -s 7
