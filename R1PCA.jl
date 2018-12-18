include(homedir()*"/Gits/PCA/cPCA.jl");

function R1PCA(inpX::Array{Float64,2}; U₀=nothing, k=Inf)
  function wH(x::Array{Float64,1},c::Float64, U::Array{Float64,2})
    norm(x-U*U'*x)<=c ? 1 : c/norm(x-U*U'*x)
  end
  function wH(x::Array{Float64,1},c::Float64, U::Array{Float64,1})
    norm(x-U*U'*x)<=c ? 1 : c/norm(x-U*U'*x)
  end
  const Loop = 100;
  const ϵ = 1e-8;

  (N,p)=size(inpX);
  meanv=mean(inpX,1)'|>vec;
  X=inpX.-meanv';
  if k==Inf
    k=p;
  end
  if U₀ == nothing
    U₀=cPCA(X,iscentered=true)[2][:,p-k+1:p];
  end

  s=map(i->sqrt(max(dot(X[i,:],(eye(p)-U₀*U₀')*X[i,:]),0)), 1:N);
  c=median(s);

  U=U₀;
  for l in 1:Loop
    q=map(i->wH(X[i,:],c,U),1:N);
    W=X'*(X.*q);
    Uˢ=qr(W*U)[1];
    if norm(Uˢ-U)<ϵ
      U=Uˢ;
      break;
    else
      U=Uˢ;
    end
  end

  U
end
