for r in run00*;
do pushd $r;
	./runme.sh;
	popd;
done
