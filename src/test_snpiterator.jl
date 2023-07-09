test_snparray_iterator(d)
d_iterator = SnpArrayIterator(d)
@test length(d_iterator) == 6
@test iterate(d_iterator,7) === nothing
@test_throws BoundsError iterate(d_iterator,0)

f = create_dummy("test2", 4, 7)
@show f
test_snparray_iterator(f)
f_iterator = SnpArrayIterator(f)
@test length(f_iterator) == 6
@test iterate(f_iterator,5) == (["1", "1", "1", "1", "1", "1", "1"],6)
@test iterate(f_iterator,6) == (["2", "2", "2", "2", "2", "2", "2"],7)
@test_throws BoundsError iterate(f_iterator,0)

g = create_dummy("test3",1,1)
@show g
test_snparray_iterator(g)
g_iterator = SnpArrayIterator(g)
@test iterate(g_iterator,1) == (["1"],2)
@test_throws BoundsError iterate(g_iterator,0)

