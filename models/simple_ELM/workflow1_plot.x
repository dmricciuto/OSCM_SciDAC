#!/bin/bash -e

i=1
while read  sitename; do
    echo $sitename
    transpose_file.x sitecoords_select.dat | awk '{print $i}' i=$i | transpose_file.x > tmppp
    lon=`awk '{print $1}' tmppp`
    lat=`awk '{print $2}' tmppp`
    echo $lon, $lat
    /Users/ksargsy/research/merr_dist/run/merr_main/pick_ind.py $lon $lat 0 None 0 -x xdata.txt -i ind_plot_${sitename}.dat
    if [[ $sitename == "US-PFa" ]]; then
        head -n12 ind_plot_${sitename}.dat > ind_plot.dat
    elif [[ $sitename == "US-WCr" ]]; then
        tail -n12 ind_plot_${sitename}.dat > ind_plot.dat
    else
        cp ind_plot_${sitename}.dat ind_plot.dat
    fi

    /Users/ksargsy/research/merr_dist/run/merr_main/plot_shade.py -i ind_plot.dat -f fixindnom.dat -d ydata.txt  -e 10 -y ytrain.txt -p pchain.dat -t $sitename
    #prep xcond_names
    /Users/ksargsy/research/merr_dist/run/merr_main/plot_fit1.py -i ind_plot.dat -x xdata.txt -y ydata.txt -m ytrain.txt -s surr_errors.dat -f fixindnom.dat -p pchain.dat -c xcond_names.txt
    cp fit.eps fit_${sitename}_red.eps
    /Users/ksargsy/research/merr_dist/run/merr_main/plot_fit1.py -i ind_plot.dat -x xdata.txt -y ydata.txt -s surr_errors.dat -f fixindnom.dat -p pchain.dat -c xcond_names.txt -t $sitename
    cp fit.eps fit_${sitename}.eps

    i=`expr $i + 1`
done < sitenames_select.dat



#/Users/ksargsy/research/merr_dist/run/merr_main/pick_ind.py -82.25 48.25 0 None 0 -x xdata_all.txt -i ind_plot.dat



