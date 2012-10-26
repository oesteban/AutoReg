#!/bin/sh
#

if [ $# -lt 3 ]
then
   echo "Monica's TFM registration script. Error: Bad arguments\n"
   echo "Usage: $0 id data_folder output_prefix"
   exit 1
fi

#### Escogemos las segmentaciones a registrar. Imagen fija: na05.nii
id=$1
data_folder=${2%/}
output_prefix=$3

mkdir -p `dirname $output_prefix`


FIX=${data_folder}/na05.nii
MOV=${data_folder}/na${id}.nii
SEG=${data_folder}/na${id}_seg.nii

OUTPUTNAME="${output_prefix}reg${id}-5MI"
OUTPUT="-o ${OUTPUTNAME}"
DIM="3" #Dimensión de la imagen


#Diversas iteraciones--> 3 niveles de optimización
IT1="-i 30x20x10" 
IT2="-i 50x30x10"

#Distintos métodos de transformación
TSYN1="-t Syn[1.0] -r Gauss[6.25,0.25]"
TSYN2="-t SyN[1.5] -r Gauss[3,0.5]"

#Diversas métricas 
METCC="-m CC[${FIX},${MOV},1,4]"
METPR="-m PR[${FIX},${MOV},1,2]"
METMI="-m MI[${FIX},${MOV},1,32]"


# Definición Parámetros WARP
#INVW=" -i ${OUTPUTNAME}Affine.txt ${OUTPUTNAME}InverseWarp.nii.gz "
FWDW="${OUTPUTNAME}Warp.nii.gz ${OUTPUTNAME}Affine.txt"

# Elegimos el grado de libertad,la transformación y la métrica deseada
ITS=$IT2
TRAN=$TSYN2
MET=$METMI

echo "------------------------------------------------------------------------------"
echo "                                PARAMETROS                                    "
echo "------------------------------------------------------------------------------"
echo " EJECUTAMOS ANTS CON LOS SIGUIENTES PARAMETROS"
echo "....SALIDA:"
eval echo $OUTPUT  
echo "....ITERACIONES":
eval echo $ITS
echo "....TRANSFORMACION:"
eval echo $TRAN
echo "....METRICA:"
eval echo $MET 
echo "------------------------------------------------------------------------------"
echo "                               FIN PARAMETROS ANTS                            "
echo "------------------------------------------------------------------------------" 


# Ejecutamos el registro ANTS
# Example usage ANTS ImageDimension -m MI[fixedimage.nii.gz,movingimage.nii.gz,1,32] -o Outputfname.nii.gz -i 30x20x0 -r Gauss[3,1] -t Elast[3]
ants_cmd="/usr/bin/time --output=${output_prefix}time_${id}.txt ANTS $DIM $MET $OUTPUT $ITS $TRAN > ${output_prefix}ants_${id}.log"

echo $ants_cmd
eval $ants_cmd
#/usr/bin/time --output=${output_prefix}time_${id}.txt "$ants_cmd"


echo "--------------------------------------------------------"
echo " FIN DE LA TRANSFORMACIÓN "
echo "--------------------------------------------------------"

# Aplicamos método WARP a la transformación resultante
# cd /usr/local/bin
# Forward Warp
fixed_warpcmd="WarpImageMultiTransform ${DIM} ${MOV} ${OUTPUTNAME}na${id}to05.nii.gz -R ${FIX} $FWDW"
echo $fixed_warpcmd
eval $fixed_warpcmd

# Aplicamos el Warp resultante a la segmentación
# Ejemplo: WarpImageMultiTransform 3 ../data/na01_seg.nii na01_segtoFIXtrilineal.nii.gz -R ../data/na05 outWarp.nii.gz outAffine.txt
seg_warpcmd="WarpImageMultiTransform ${DIM} ${SEG} ${OUTPUTNAME}seg${id}to05.nii.gz -R ${FIX} $FWDW --use-NN"
echo $seg_warpcmd
eval $seg_warpcmd


#AÑADIR NN nearest neighbor, pq sino esta interpolando y no molaaaa..histograma raroo!

echo "--------------------------------------------------------"
echo " PROCESO DE REGISTRO COMPLETADO "

echo "--------------------------------------------------------"

echo "--------Eliminando imagenes invertidas--------------------"

#rm -rf ${OUTPUTNAME}InverseWarp.nii.gz

echo "IMAGEN REGISTRADA: ${OUTPUTNAME}MOVtoFIX.nii.gz"
echo "IMAGEN SEGMENTADA FINAL REGISTRADA: ${OUTPUTNAME}SEG${id}_FIX.nii.gz"



# Guardamos los ficheros generados en una carpeta con el mismo nombre
# mkdir /home/mferrer/Escritorio/Registros/${OUTPUTNAME}
# mv /usr/lib/ants/${OUTPUTNAME}MOVtoFIX.nii.gz ${OUTPUTNAME}Warp.nii.gz ${OUTPUTNAME}Affine.txt /home/mferrer/Escritorio/Registros/${OUTPUTNAME}
