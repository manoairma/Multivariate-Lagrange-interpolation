\\ Génère les vecteurs dont la somme est inferieur ou egal à n
\\ Chaque ei varie de 0 à n
gen_exponents(p,n) = {
  my(L = List());
  \\ forvec itère sur tous les vecteurs [e1, ..., ep] où chaque composante ei est entre 0 et n
  forvec(v = vector(p, i, [0, n]),
    if(vecsum(v) <= n, listput(L, v))
  );
  return(Vec(L));
}
\\Interpolation de Lagrange d'une fonction a p variable par un polynome de degree total n
\\P:liste des points 
\\F:valeur de la fonction
\\X:variable
\\n: degre de l'interpolant

lagrange_multi(P,F,X,n)={
  my(p=#X);
  my(E=gen_exponents(p,n));
  my(k=#E);

\\verification du nombre de points
  if(k!=#P||k!=#F,
    error("le nombre de points doit etre binomial(p+n,n)")
  );
  
\\calcul de la matrice
  my(M=matrix(k,k));
  for(i=1,k,
    for(j=1,k,
        M[i,k-j+1]=prod(v=1,p,P[i][v]^E[j][v]);
    );
  );
  
  \\calcul de déterminant de M
  my(Delta= matdet(M));
    if(delta==0,
    error("il y a des points qui se répètent");
 );
 
 \\ Construction du vecteur x qui va remplacer la iéme ligne de M
  my(x= vector(k, j, prod(v=1, p, X[v]^E[k-j+1][v])));

  \\ Calcul de P(X) = Somme Fi * Li(X)  avec li(X) = det(Mi)/det(M)
  my(P_X = 0);
  for(i=1, k,
    my(Mi = M);
    Mi[i, ] = x; \\ On remplace la ligne i par la variable standard x
    P_x += F[i] * matdet(Mi) / Delta;
  );
 
  return([M , Delta , P_x]);

}

lagrange_multi([[0,0],[0,1],[1,1]],[1,2,3],[x_1,x_2],1);

