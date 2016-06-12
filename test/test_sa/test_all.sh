for func in load summarize grm mom pca pca_sp
do
    julia --precompiled=yes test_"$func".jl > ./result/test_"$func"_result.txt
    sleep 60
done
