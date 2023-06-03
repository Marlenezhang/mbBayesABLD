#!/usr/bin/bash

## 对提供的任意多个群体进行PCA分析并作图
## 需要用到的软件：R、Plink1.9
## 需要提供的信息，各个群体的plink文件前缀（plink文件需转为fam、bim、bed的二进制格式）
## 存在不在当前工作路径的文件时，提供全路径，但不用带后缀(.fam)，如/home/user/pop1

## 用法：./pca_multi_pop.sh -P "prefix1 prefix2 ..." -N "N1 N2 ..."
# -P 需要作图的群体的plink文件路径（需为二进制文件，且路径中文件名不用加.fam .bim .bed等后缀）
# -N 每个群体在PCA作图中的群体标识符 (两个参数均需加双引号)
# -O 输出文件名
# -h 帮助

###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long pre_list:,fids:,out:,help,fid,plot \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
          --pre_list )   pre_list="$2";  shift 2 ;;  ## plink文件前缀，如"/public/home/popA /public/home/popB"
          --fids )       fids="$2"; 	   shift 2 ;;	 ## 品种（群体）标识符，如"popA popB"
          --out )        out="$2" ;      shift 2 ;;  ## 输出文件名
          --fid )        fid=true;       shift   ;;  ## bfile中提供了品种（群体）标识符(即fid)
          --plot )       plot=true;      shift   ;;  ## 根据PCA结果作图
    -h | --help )        grep " shift " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

# ## 加载需要的软件
# module load R/4.1.0
# module load PLINK/1.90

## 脚本
code=${code:=/BIGDATA2/cau_jfliu_2/liwn/code}      ## 脚本存放目录
pca_plot=${code}/R/plot/PCA_plot.R

## 检查必要参数
[[ ! ${pre_list} ]] && echo "Open script $0 to view instructions" && exit 1

## 只提供了一个文件，则根据fam文件中的fid进行分组
files_num=$(echo ${pre_list} | tr " " "\n" | wc -l)
if [[ ${files_num} -le 1 ]]; then
  unset fids
  if [[ ! ${fid} ]]; then
    echo "please provide the --fid option if the .ped(.fam) file contains fid"
    exit 1
  fi

  ## 检查文件是否存在
  if [[ ! -f ${pre_list}.fam ]]; then
    if [[ -f ${pre_list}.ped ]]; then
      ## 格式转换
      plink --file ${pre_list} --make-bed --out ${pre_list} > plink.log
    else
      echo "plink file ${pre_list} not find, please check! "
      exit 1
    fi
  fi

  fid_uniq=$(awk '{print $1}' ${pre_list}.fam | sort | uniq)

  ## 提取不同品种的基因型信息
  for name in ${fid_uniq}; do
    awk ${name} > ${name}_fid.txt
    plink --bfile ${pre_list} --keep-fam ${name}_fid.txt --chr-set 30 --make-bed --out ${name}
  done

  pre_list=("${fid_uniq}")
  # Name_list="${pre_list[*]}"
fi

## 获取每个群体的标识符（family id）
# pre_list=("${pre_list[@]}")
IFS=' ' read -ra pre_list <<< "${pre_list[@]}"
index=0
for prefix in "${pre_list[@]}"; do # prefix=${pre_list[0]}
  ((index++))

  ## 检查文件是否存在
  if [[ ! -f ${prefix}.fam ]]; then
    if [[ -f ${prefix}.ped ]]; then
      ## 格式转换
      plink --file ${prefix} --make-bed --out ${prefix} > plink.log
    else
      echo "plink file ${prefix} not find, please check! "
      exit 1
    fi
  fi

  ## 群体id
  if [[ ${fids[*]} ]]; then
    fidi=$(echo "${fids[*]}" | cut -d ' ' -f ${index})
  else
    ## 提取第一行iid和fid，判断fid是否用iid表示或全为0，是的话则重新给个fid（如pop1）
    iid=$(head -n 1 ${prefix}.fam | awk '{print $2}')
    fidi=$(head -n 1 ${prefix}.fam | awk '{print $1}')
    if [[ ${fidi} == "${iid}" || ${fidi} == '0' ]]; then
      fidi=pop${index}
    fi
    # Name_list="${Name_list} ${fidi}"
  fi

  ## 输出iid和fid匹配表
  if [[ ${index} -eq 1 ]]; then
    awk '{print $2,"'${fidi}'"}' ${prefix}.fam > iid_fid_pca.txt
  else
    awk '{print $2,"'${fidi}'"}' ${prefix}.fam >> iid_fid_pca.txt

    ## 合并群体需要的文件名列表
    if [[ ${index} -ge 2 ]]; then
      echo ${prefix}.bed ${prefix}.bim ${prefix}.fam >> merge_list.txt
    else
      [[ -f merge_list.txt ]] && rm merge_list.txt
    fi
  fi
done

## 输出文件名
[[ ! ${out} ]] && out=$(echo "${pre_list[@]}" | tr " " "_")_pca.txt

## 合并所有群体plink文件
prefix1=$(echo "${pre_list[@]}" | cut -d ' ' -f 1)
num=$(echo "${pre_list[@]}" | wc -w)
plink --bfile ${prefix1} --merge-list merge_list.txt --chr-set 30 --make-bed --out pop${num}_merge > plink.log

## 存在多等位基因位点
if [[ $? -ne 0 ]]; then
  echo " remove all offending variants"
  plink --bfile source1 --exclude merged.missnp --make-bed --out source1_tmp > plink.log
  plink --bfile source2 --exclude merged.missnp --make-bed --out source2_tmp
  plink --bfile source1_tmp --bmerge source2_tmp --make-bed --out merged
  rm source1_tmp.*
  rm source2_tmp.*
fi

## pca计算
plink --bfile pop${num}_merge --chr-set 30 --pca 3 header --out pop${num}_merge > plink.log

## 匹配群体标识符
if [[ ! ${fid} ]]; then
  sort -k 1n iid_fid_pca.txt -o iid_fid_pca.txt2
  sed '1d' pop${num}_merge.eigenvec | sort -k 2n > pop${num}_merge.eigenvec2
  join -1 2 -2 1 pop${num}_merge.eigenvec2 iid_fid_pca.txt2 > "${out}"
else
  cp pop${num}_merge.eigenvec "${out}" 
fi

## 作图
if [[ ${plot} ]]; then
  $pca_plot --eigv ${out}
fi

## 删除中间文件
rm pop${num}_merge.*
rm merge_list.*
rm iid_fid_pca.*

## debug
# cd /BIGDATA2/cau_jfliu_2/liwn/mbGS/ChipSim/Lin2022
pre_list="/BIGDATA2/cau_jfliu_2/liwn/mbGS/Real/Lee2019/Lee2019q"
pre_list="STHqm SNTqm QBqm"
# Name_list="A B C D"
# out="breedA_breedB_pca.txt"
fid=true
plot=true