
for index in {1..5..1}
do
	echo $index
	sed -e "s/XXX/$index/g" < rmsd.config > temp.config
	python3 rmsd.ref.py temp.config > model_$index.out &
	sleep 10s
done

wait

python3 rmsd_plotting.py model_1.00.19.hgf_wt.rmsd selection_list.out model_1 &
python3 rmsd_plotting.py model_2.00.19.hgf_wt.rmsd selection_list.out model_2 &
python3 rmsd_plotting.py model_3.00.19.hgf_wt.rmsd selection_list.out model_3 &
python3 rmsd_plotting.py model_4.00.19.hgf_wt.rmsd selection_list.out model_4 &
python3 rmsd_plotting.py model_5.00.19.hgf_wt.rmsd selection_list.out model_5 &

wait

