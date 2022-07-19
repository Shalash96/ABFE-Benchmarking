extension="dat" # the extension of the files you want to search in
word="lig"      # the word that you want to search for
for file in $ls *.$extension
do
    grep -i "lig" $file > $(basename $file .$extension)_lig.$extension
    if [ ! -s $(basename $file .$extension)_${word}.$extension ]
    then
        rm $(basename $file .$extension)_${word}.$extension
    fi
done

