using LinearAlgebra, Statistics

function R1PCA(inpX, k; centered=2)
  δ=1e-10;
  ϵ=1e-10;
  Loop=100;
  function GM(X)
    (N,d)=size(X);
    cⁿ=zeros(d);
    cˢ=zeros(d);

    cⁿ.=mean(X,dims=1)|>vec;
    for l in 1:Loop
      ds=map(i->1/max(δ, norm(X[i,:]-cⁿ)), 1:N);
      cˢ.=(mean(X .* ds, dims=1) |>vec) / sum(ds);

      if(norm(cˢ-cⁿ) < ϵ)
        cⁿ.=cˢ;
        break;
      else
        cⁿ=cˢ;
      end
    end
    return(cⁿ);
  end
  function wH(x, U, c)
    norm(x-U*U'*x)<=c ? 1 : c/norm(x-U*U'*x)
  end

  (N,d)=size(inpX);
  X=zeros(N,d);

  if centered==1
    μ=GM(inpX);
    X.=inpX .- μ';
  elseif centered==2
    μ=mean(inpX,dims=1)|>vec;
    X.=inpX .- μ';
  end

  (λs, U₀) = eigen(X'*X);
  U₀.=U₀[:,sortperm(λs, rev=true)];
  U=zeros(d,k);
  U.=U₀[:,1:k];

  s=map(
    i->sqrt(
      max(
        dot(X[i,:], (Matrix(1.0I, d, d)-U*U')*X[i,:]),
        0
      )
    ),
    1:N
  );
  c=median(s);

  for l in 1:Loop
    q=map(i->wH(X[i,:],U, c),1:N);
    W=X'*(X.*q);
    (Uˢ,R)=qr(W*U);
    if norm(Uˢ-U) < ϵ
      U=Uˢ;
      break;
    else
      U=Uˢ;
    end
  end

  U
end
