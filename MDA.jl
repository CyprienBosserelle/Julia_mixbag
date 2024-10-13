

function normalise(x; minx=NaN,maxx=NaN)

	if isnan(minx)
		minx=minimum(x);
	end

	if isnan(maxx)
		maxx=maximum(x);
	end

	norm = (x .- minx) ./ (maxx - minx)



	return norm
end


"""
### Maximum dissimilarity algorithm MDA

here we assume data is m x n matrix where m (rows) is the different dimensions and n (columns) are the different scattered data.
seed is the initial seed index for MDA (needs to be <= n)
"""
function mda(data::Matrix{Float64},ncentres;isdeg=(), seed=1)

	# Normalise data
	m,n=size(data);
	morndata=zeros(size(data));

	centerindex=Vector{Int64}(undef, 0)

	for i=1:m
		if any(i .== isdeg)
			morndata[i,:] = data[i,:] .* pi ./ 180.0
		else
			morndata[i,:]=normalise(data[i,:])
		end
	end

	dataindex=collect(1:n);

	push!(centerindex,seed)

	popat!(dataindex,seed);

	indxtokeep=1:n.!=seed


	lastcenter=morndata[:,seed];

	#println(lastcenter)

	subset=morndata[:,indxtokeep]
	lastdist=colwise(SqEuclidean(), lastcenter, subset)


	while length(centerindex)<ncentres
		r = colwise(SqEuclidean(), lastcenter, subset)
		lastdist=min.(r,lastdist);

		idx=argmax(lastdist);

		#println(r)
		push!(centerindex,dataindex[idx]);

		popat!(dataindex,idx)

		popat!(lastdist,idx)

		p,q=size(subset)

		indxtokeep=1:q.!=idx

		lastcenter=subset[:,idx];

		subset=subset[:,indxtokeep]
	end

	return centerindex

end
