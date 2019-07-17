using LinearAlgebra, Statistics
include(homedir()*"/Gits/sources/GeometricMedian.jl");

function R1PCA(inpX::Array{Float64,2}, k::Int64; isCentered=false)
  δ=1e-10;
  ϵ=1e-10;
  Loop=100;
  function wH(x, U, c)
    norm(x-U*U'*x)<=c ? 1 : c/norm(x-U*U'*x)
  end

  (N,d)=size(inpX);
  X=zeros(N,d);
  gm=zeros(d);

  if !isCentered
    gm=GM(inpX);
  end
  X.=inpX.-gm';

  eigret = X'*X |> Symmetric |> eigen;
  λs = eigret.values;
  U₀ = eigret.vectors[:, sortperm(λs, rev=true)];
  U=zeros(d,k);
  U.=U₀[:,1:k];

  s=map(
    i->sqrt(
      max(0, dot(X[i,:], (Matrix(1.0I, d, d)-U*U')*X[i,:]))
    ),1:N
  );
  c=median(s);

  for l in 1:Loop
    q=map(i->wH(X[i,:],U, c),1:N);
    W=X'*(X.*q);
    Uˢ=qr(W*U).Q |> Matrix;
    if norm(Uˢ-U) < ϵ
      U=Uˢ;
      break;
    else
      U=Uˢ;
    end
  end

  return (gm, U);
end
