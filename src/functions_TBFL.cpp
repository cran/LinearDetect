// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

//' soft threshold function
//'
//' @param L input matrix
//' @param lambda threshold parameter  
//' @return thresholded matrix L
// [[Rcpp::export]]
arma::mat soft_full(arma::mat L, double lambda){
	int height = L.n_rows;
	int width = L.n_cols;
 	for(int y = 0; y < height; y++){
 		for(int x = 0; x < width; x++) {
            if (L(y, x) > lambda)
            	L(y, x) = L(y, x) - lambda;
            else if (L(y, x) <  - lambda)
                L(y, x) = L(y, x) + lambda;
            else
                L(y, x) = 0.0;
        }

 	} 
 	return L;
}

//[[Rcpp::export]]
List lm_break_fit_block(NumericMatrix data_y, NumericMatrix data_x, double lambda, double lambda2, int max_iteration, double tol , NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index ){

    
    int p_y = data_y.ncol(); int p_x = data_x.ncol(); int T = data_y.nrow(); 
    int n_new = blocks.size() - 1;
    int n = T;


    arma::mat data_y_m(data_y.begin(), T, p_y);
    arma::mat data_x_m(data_x.begin(), T, p_x);
    List y_b(n_new);

    for(int i = 0; i < n_new; i++) {
    	y_b[i] = data_y_m( span(blocks[i]-1,  blocks[i+1]-2 ), span::all) ;
    }

    
    arma::mat X( p_x, n);
    for(int i = 0; i < T; i++) {
            X(span::all, i )= data_x_m( i , span::all ).t();
    }

    List X_b(n_new);   
    for(int i = 0; i < n_new; i++) {
        X_b[i] = X(span::all, span(blocks[i]-1, blocks[i+1]-2 ));
    }

    int cv_l = cv_index.size();
    if( cv_l >0){
    	for(int t_1 = 0; t_1 < cv_l; t_1 ++){
    		arma::mat yb_temp = y_b[cv_index[t_1]-1]; 
    		arma::mat Xb_temp = X_b[cv_index[t_1]-1]; 
    		int tt =  yb_temp.n_rows;
            if( tt > 1 ){
                y_b[cv_index[t_1]-1] = yb_temp(span(0, tt-2), span::all );
                X_b[cv_index[t_1]-1] = Xb_temp( span::all, span(0, tt-2));
            }
            else{
                Rcout << "Wrong!";
            }          

    	}
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
    	arma::mat yb_temp = y_b[j];
    	arma::mat Xb_temp = X_b[j]; 
    	C[j] = Xb_temp * yb_temp;
    }  
    
    arma::mat C_sum (p_x*n_new, p_y, fill::zeros);
    //Rcout << C_sum(1,1);
    arma::mat C_temp = C[0];
    C_sum(span(0, p_x - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
    	arma::mat C_temp = C[i-1];
    	//Rcout <<  size(C_temp);
    	C_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = C_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + C_temp; 
    }

    arma::mat C_sum_new (p_x*n_new, p_y, fill::zeros);
    C_sum_new( span(   0, p_x - 1 ), span::all) =  C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
    	C_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - C_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
    }


    List D(n_new);
    for(int j = 0; j < n_new; j++){
    	mat Xb_temp = X_b[j]; 
    	D[j] = Xb_temp * Xb_temp.t();
    }

    arma::mat D_sum (p_x*n_new, p_x, fill::zeros);
    arma::mat D_temp = D[0];
    D_sum(span(0, p_x - 1) , span::all) = D_temp;
    for(int i = 2; i <= n_new ; i++){
    	arma::mat D_temp = D[i-1];
    	//Rcout <<  size(C_temp);
    	D_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = D_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + D_temp; 
    }
    arma::mat D_sum_new (p_x*n_new, p_x, fill::zeros);
    D_sum_new( span(  0, p_x - 1 ), span::all) =  D_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
    	D_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = D_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - D_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
    }

    // Rcout << D_sum_new(  span(0, 9), span(0, 9) ) << "\n";

    arma::mat D_sum_new_inv (p_x*n_new, p_x, fill::zeros);

    for(int i = 1; i <= n_new ; i++){
    	// Rcout << D_sum_new ( span((i-1)*p_x, (i*p_x)-1), span::all);
        vec eigval =  eig_sym(D_sum_new ( span((i-1)*p_x, (i*p_x)-1), span::all)  );
        // if(i == 1){
        //     Rcout << "eigen values are" << std::endl << eigval << std::endl;
        // }
        double add_pd = 0.0;
        double min_eigen = eigval(0);
        if(min_eigen <= 0){
            Rprintf("not invertiable! adding noise!");
            if(min_eigen < 0 ){
                // Rcout << i << "\n";
                // Rcout << fabs( min_eigen) << "\n"; 
            	// add_pd = fabs( min_eigen) + pow(10, -3);
            	add_pd = fabs( min_eigen) + pow(10.0, -6.0);
                // Rcout << fabs( add_pd) << "\n"; 
                // add_pd = pow(10, -3); 
            }
            else
            	add_pd = pow(10.0, -6.0);           
        }        
        arma::mat noise(p_x, p_x, fill::eye);
        D_sum_new_inv(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = inv( D_sum_new( span((i-1)*p_x, (i*p_x) -1  ), span::all  ) +  add_pd*noise );
    }

    // Rcout << D_sum_new_inv.n_rows << "\n";
    // Rcout << D_sum_new_inv.n_cols << "\n";

    // Rcout << "D sum inv matrix is" << std::endl << D_sum_new_inv( span ((1-1)*p_x , 1*p_x-1), span::all ) << std::endl;


    arma::mat phi_hat(initial_phi.begin(), p_y, p_x*n_new);
    arma::mat phi_new (p_y, p_x*n_new, fill::zeros);

    int flag = 0;
    int l = 2;
    
    while( l < max_iteration){
    	if(  l == floor(0.5 * max_iteration)  ){
    		tol = (2)*tol;
    	}
    	if( l  == floor( 0.75 * max_iteration )) {
    		tol = (4/2)*tol;
    	}
        l = l+1; 
        arma::mat phi_compare = phi_hat;

        for(int i = 1; i <= n_new ; i++){
        	List E(n_new);
        	for(int j = 1; j <= n_new ; j++){       		
        		E[j-1] = D_sum_new( span(  ( std::max( j,i )-1)*p_x  , std::max(j,i )*p_x  -1 ), span::all)  * phi_hat(span::all , span( (j-1)*p_x, j*p_x -1  )).t() ;    			
        	}

        	
        	arma::mat E_new = E[0];
    		for ( int g = 1; g < n_new; g++ ) {
    			arma::mat E_temp = E[g];
    			E_new = E_new + E_temp;
    		}
    		

    		E_new =  E_new - D_sum_new(  span( (i-1)*p_x,  i*p_x -1  ), span::all)  * phi_hat( span::all, span( (i-1)*p_x, i*p_x -1) ).t();
            // Rcout << E_new;

            // Rcout <<  C_sum_new(  span(  (1-1)*p_x , 1*p_x -1 ), span::all );
    		arma::mat S =  C_sum_new(  span(  (i-1)*p_x , i*p_x -1 ), span::all )  - E_new;

            // if(i == 1){
            //     Rcout << "S matrix was" << std::endl << S << std::endl;

            // }

            
   		

    		S = soft_full(S, lambda);


            // if(i == 1){
            //     Rcout << "S matrix is" << std::endl << S << std::endl;
            // }


    		arma::mat phi_temp = D_sum_new_inv( span ((i-1)*p_x , i*p_x-1), span::all )  *  S;

            // if(i == 1){
            //     Rcout << "D sum inv matrix is" << std::endl << D_sum_new_inv( span (0, 3), span (240, 249)) << std::endl;
            // }
    		
    		phi_temp = phi_temp.t();


    		phi_hat( span::all ,   span( (i-1)*p_x, i*p_x -1)  ) = phi_temp;
            phi_new( span::all ,   span( (i-1)*p_x, i*p_x -1)  ) = phi_temp;
                                                                  

        }

        // Rcout << phi_hat;

        arma::mat phi_temp_soft(p_y,p_x*n_new, fill::zeros);
        phi_temp_soft( span::all , span(0, p_x-1) ) = soft_full(phi_hat( span::all ,  span(0, p_x-1) ), lambda2);
        arma::mat temp_1 = phi_hat( span::all ,  span(0, p_x-1) );
        for(int z_1 = 2; z_1 <= n_new; z_1 ++){
            arma::mat temp_2 = temp_1 + phi_hat( span::all ,  span((z_1-1)*p_x, z_1*p_x-1) );
            arma::mat temp_1_soft = soft_full(temp_1,lambda2);
            arma::mat temp_2_soft = soft_full(temp_2,lambda2);
            phi_temp_soft(   span::all ,  span( (z_1-1)*p_x , z_1*p_x  -1 ) )= temp_2_soft - temp_1_soft;
            temp_1 = temp_2;
        }
        phi_new  = phi_temp_soft;



        arma::mat abs_temp = abs(phi_new - phi_compare);
        double max_temp = abs_temp.max();
        if ( max_temp < tol) {
            break;
        } 
        if ( max_temp > tol ) {
        	phi_hat = phi_new; 
            // Rprintf( "%f \n", max_temp);

        }
        if (  max_temp > pow(10.0, 5.0)) {
            flag = 1;
            break;
        }

                
	}
	

    return List::create(Named("phi.hat")= phi_hat, Named("flag")= flag);	
}



//[[Rcpp::export]]
arma::mat lm_break_fit_block_new(arma::mat data_y, arma::mat data_x, double lambda, double lambda2, int max_iteration, double tol , NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index ){
    
    int p_y = data_y.n_cols; int p_x = data_x.n_cols; int T = data_y.n_rows; 
    int n_new = blocks.size() - 1;
    int n = T;


    arma::mat data_y_m(data_y.begin(), T, p_y);
    arma::mat data_x_m(data_x.begin(), T, p_x);
    List y_b(n_new);

    for(int i = 0; i < n_new; i++) {
    	y_b[i] = data_y_m( span(blocks[i]-1,  blocks[i+1]-2 ), span::all) ;
    }

    
    arma::mat X( p_x, n);
    for(int i = 0; i < T; i++) {
            X(span::all, i )= data_x_m( i , span::all ).t();
    }

    List X_b(n_new);   
    for(int i = 0; i < n_new; i++) {
        X_b[i] = X(span::all, span(blocks[i]-1, blocks[i+1]-2 ));
    }

    int cv_l = cv_index.size();
    if( cv_l >0){
    	for(int t_1 = 0; t_1 < cv_l; t_1 ++){
    		arma::mat yb_temp = y_b[cv_index[t_1]-1]; 
    		arma::mat Xb_temp = X_b[cv_index[t_1]-1]; 
    		int tt =  yb_temp.n_rows;
            if( tt > 1 ){
                y_b[cv_index[t_1]-1] = yb_temp(span(0, tt-2), span::all );
                X_b[cv_index[t_1]-1] = Xb_temp( span::all, span(0, tt-2));
            }
            else{
                Rcout << "Wrong!";
            }          

    	}
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
    	arma::mat yb_temp = y_b[j];
    	arma::mat Xb_temp = X_b[j]; 
    	C[j] = Xb_temp * yb_temp;
    }  
    
    arma::mat C_sum (p_x*n_new, p_y, fill::zeros);
    //Rcout << C_sum(1,1);
    arma::mat C_temp = C[0];
    C_sum(span(0, p_x - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
    	mat C_temp = C[i-1];
    	//Rcout <<  size(C_temp);
    	C_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = C_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + C_temp; 
    }

    arma::mat C_sum_new (p_x*n_new, p_y, fill::zeros);
    C_sum_new( span(   0, p_x - 1 ), span::all) =  C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
    	C_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - C_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
    }


    List D(n_new);
    for(int j = 0; j < n_new; j++){
    	arma::mat Xb_temp = X_b[j]; 
    	D[j] = Xb_temp * Xb_temp.t();
    }

    arma::mat D_sum (p_x*n_new, p_x, fill::zeros);
    arma::mat D_temp = D[0];
    D_sum(span(0, p_x - 1) , span::all) = D_temp;
    for(int i = 2; i <= n_new ; i++){
    	arma::mat D_temp = D[i-1];
    	//Rcout <<  size(C_temp);
    	D_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = D_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + D_temp; 
    }
    arma::mat D_sum_new (p_x*n_new, p_x, fill::zeros);
    D_sum_new( span(  0, p_x - 1 ), span::all) =  D_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
    	D_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = D_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - D_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
    }

    // Rcout << D_sum_new(  span(0, 9), span(0, 9) ) << "\n";

    arma::mat D_sum_new_inv (p_x*n_new, p_x, fill::zeros);

    for(int i = 1; i <= n_new ; i++){
    	// Rcout << D_sum_new ( span((i-1)*p_x, (i*p_x)-1), span::all);
        vec eigval =  eig_sym(D_sum_new ( span((i-1)*p_x, (i*p_x)-1), span::all)  );
        // if(i == 1){
        //     Rcout << "eigen values are" << std::endl << eigval << std::endl;
        // }
        double add_pd = 0.0;
        double min_eigen = eigval(0);
        if(min_eigen <= 0){
            Rprintf("not invertiable! adding noise!");
            // Rcout << i << "\n"; 
            // Rcout << fabs( min_eigen) << "\n";
            if(min_eigen < 0 ){
            	add_pd = (10)*fabs( min_eigen);
            }
            else
            	add_pd = pow(10.0, -6.0);          
        }    
        arma::mat noise(p_x, p_x, fill::eye);
        D_sum_new_inv(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = inv( D_sum_new( span((i-1)*p_x, (i*p_x) -1  ), span::all  ) +  add_pd*noise );
    }


    // Rcout << "D sum inv matrix is" << std::endl << D_sum_new_inv( span ((1-1)*p_x , 1*p_x-1), span::all ) << std::endl;


    arma::mat phi_hat(initial_phi.begin(), p_y, p_x*n_new);
    arma::mat phi_new (p_y, p_x*n_new, fill::zeros);

    int flag = 0;
    int l = 2;
    
    while( l < max_iteration){
    	if(  l == floor(0.5 * max_iteration)  ){
    		tol = (2)*tol;
    	}
    	if( l  == floor( 0.75 * max_iteration )) {
    		tol = (4/2)*tol;
    	}
        l = l+1; 
        arma::mat phi_compare = phi_hat;

        for(int i = 1; i <= n_new ; i++){
        	List E(n_new);
        	for(int j = 1; j <= n_new ; j++){       		
        		E[j-1] = D_sum_new( span(  ( std::max( j,i )-1)*p_x  , std::max(j,i )*p_x  -1 ), span::all)  * phi_hat(span::all , span( (j-1)*p_x, j*p_x -1  )).t() ;    			
        	}

        	
        	arma::mat E_new = E[0];
    		for ( int g = 1; g < n_new; g++ ) {
    			arma::mat E_temp = E[g];
    			E_new = E_new + E_temp;
    		}
    		

    		E_new =  E_new - D_sum_new(  span( (i-1)*p_x,  i*p_x -1  ), span::all)  * phi_hat( span::all, span( (i-1)*p_x, i*p_x -1) ).t();
    		arma::mat S =  C_sum_new(  span(  (i-1)*p_x , i*p_x -1 ), span::all )  - E_new;
   		

    		S = soft_full(S, lambda);



    		arma::mat phi_temp = D_sum_new_inv( span ((i-1)*p_x , i*p_x-1), span::all )  *  S;


    		
    		phi_temp = phi_temp.t();


    		phi_hat( span::all ,   span( (i-1)*p_x, i*p_x -1)  ) = phi_temp;
            phi_new( span::all ,   span( (i-1)*p_x, i*p_x -1)  ) = phi_temp;
                                                                  

        }

        // Rcout << phi_hat;

        arma::mat phi_temp_soft(p_y,p_x*n_new, fill::zeros);
        phi_temp_soft( span::all , span(0, p_x-1) ) = soft_full(phi_hat( span::all ,  span(0, p_x-1) ), lambda2);
        arma::mat temp_1 = phi_hat( span::all ,  span(0, p_x-1) );
        for(int z_1 = 2; z_1 <= n_new; z_1 ++){
            arma::mat temp_2 = temp_1 + phi_hat( span::all ,  span((z_1-1)*p_x, z_1*p_x-1) );
            arma::mat temp_1_soft = soft_full(temp_1,lambda2);
            arma::mat temp_2_soft = soft_full(temp_2,lambda2);
            phi_temp_soft(   span::all ,  span( (z_1-1)*p_x , z_1*p_x  -1 ) )= temp_2_soft - temp_1_soft;
            temp_1 = temp_2;
        }
        phi_new  = phi_temp_soft;



        arma::mat abs_temp = abs(phi_new - phi_compare);
        double max_temp = abs_temp.max();
        if ( max_temp < tol) {
            break;
        } 
        if ( max_temp > tol ) {
        	phi_hat = phi_new; 
            // Rprintf( "%f \n", max_temp);

        }
        if (  max_temp > pow(10.0, 5.0)) {
            flag = 1;
            break;
        }

                
	}
	

    return phi_hat;	
}

//[[Rcpp::export]]
List ggm_break_fit_block(NumericMatrix data_y, NumericMatrix data_x, double lambda, double lambda2, int max_iteration, double tol , NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index ){
	int p_y = data_y.ncol(); int p_x = data_x.ncol(); int T = data_y.nrow(); 
    // int n = T;
    int n_new = blocks.size() - 1;
    int p_x_temp = p_x - 1;
    // int p_y_temp = 1;

    arma::mat data_y_m(data_y.begin(), T, p_y);
    arma::mat data_x_m(data_x.begin(), T, p_x);
    arma::mat phi_hat_full_all (p_y, n_new*p_x_temp, fill::zeros);
	for(int j = 0; j < p_y ; j++){
		arma::mat data_y_temp = data_y_m( span::all ,  j );
		arma::mat data_x_temp (T, p_x-1, fill::zeros);
		for(int i = 0; i < p_y; i++){
			if( i < j ){
				data_x_temp(span::all, i) = data_x_m(span::all, i);

			}
			if( i > j ){
				data_x_temp(span::all, i-1) = data_x_m(span::all, i);

			}
		}
		arma::mat test = lm_break_fit_block_new(data_y_temp, data_x_temp, lambda, lambda2, max_iteration, tol, initial_phi, blocks, cv_index);
		arma::mat phi_hat_full = test;	
		phi_hat_full_all(j , span::all ) = phi_hat_full;

	}
	return List::create(Named("phi.hat") = phi_hat_full_all, Named("flag") = T);

}



// [[Rcpp::export]]
List lm_partial_break_fit_block(NumericMatrix data_y, NumericMatrix data_x, double lambda, double lambda2, int max_iteration, double tol , NumericMatrix initial_phi, NumericMatrix initial_phi_2, NumericVector blocks, NumericVector cv_index, arma::uvec fixed_index, arma::uvec nonfixed_index  ){

    // Environment pkg = Rcpp::Environment::namespace_env("glmnet");
    // Function f_R = pkg["cv.glmnet"];

    int p_y = data_y.ncol(); int p_x_all = data_x.ncol(); int T = data_y.nrow(); 
    int n_new = blocks.size() - 1;
    int n = T;

    int fixed_l =fixed_index.size();
    int p_x = p_x_all - fixed_l;

    arma::mat data_y_m_all(data_y.begin(), T, p_y);
    arma::mat data_x_m_all(data_x.begin(), T, p_x_all);

    arma::mat data_x_m(T, p_x);
    arma::mat data_x_m_2(T, fixed_l);
    // Rcout << fixed_index-1;
    arma::uvec row_index = regspace<arma::uvec>(0,  1,  T-1);
    // Rcout << row_index;

    data_x_m = data_x_m_all.submat(row_index, nonfixed_index -1 );
    data_x_m_2 = data_x_m_all.submat(row_index,  fixed_index -1 );

    arma::mat phi_hat(initial_phi.begin(), p_y, p_x*n_new);
    arma::mat phi_hat_2(initial_phi_2.begin(), p_y, fixed_l);

    int flag = 0;
    int l = 2;
    
    while( l < max_iteration){
    	if(  l == floor(0.5 * max_iteration)  ){
    		tol = (2)*tol;
    	}
    	if( l  == floor( 0.75 * max_iteration )) {
    		tol = (4/2)*tol;
    	}
        l = l+1; 

	    arma::mat data_y_m =  data_y_m_all - data_x_m_2* phi_hat_2.t();



	    List y_b(n_new);
	    for(int i = 0; i < n_new; i++) {
	    	y_b[i] = data_y_m( span(blocks[i]-1,  blocks[i+1]-2 ), span::all) ;
	    }

	    arma::mat X( p_x, n);
	    for(int i = 0; i < T; i++) {
	            X(span::all, i )= data_x_m( i , span::all ).t();
	    }

	    List X_b(n_new);   
	    for(int i = 0; i < n_new; i++) {
	        X_b[i] = X(span::all, span(blocks[i]-1,blocks[i+1]-2 ));
	    }


	    int cv_l = cv_index.size();
	    if( cv_l >0){
	    	for(int t_1 = 0; t_1 < cv_l; t_1 ++){
	    		mat yb_temp = y_b[cv_index[t_1]-1]; 
	    		mat Xb_temp = X_b[cv_index[t_1]-1]; 
	    		int tt =  yb_temp.n_rows;
	            if( tt>1 ){
	                y_b[cv_index[t_1]-1] = yb_temp(span(0,tt-2), span::all );
	                X_b[cv_index[t_1]-1] = Xb_temp( span::all ,span(0,tt-2));
	            }
	            else{
	                Rcout << "Wrong!";
	            }          

	    	}
	    } 

	    List C(n_new);
	    for(int j = 0; j < n_new; j++){
	    	arma::mat yb_temp = y_b[j];
	    	arma::mat Xb_temp = X_b[j]; 
	    	C[j] = Xb_temp * yb_temp;
	    }  
	    
	    arma::mat C_sum (p_x*n_new, p_y, fill::zeros);
	    arma::mat C_temp = C[0];
	    C_sum(span(0, p_x - 1) , span::all) = C_temp;
	    for(int i = 2; i <= n_new ; i++){
	    	arma::mat C_temp = C[i-1];
	    	C_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = C_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + C_temp; 
	    }

	    arma::mat C_sum_new (p_x*n_new, p_y, fill::zeros);
	    C_sum_new( span(   0, p_x - 1 ), span::all) =  C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
	    for(int i = 2; i <= n_new ; i++){
	    	C_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - C_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
	    }


	    List D(n_new);
	    for(int j = 0; j < n_new; j++){
	    	arma::mat Xb_temp = X_b[j]; 
	    	D[j] = Xb_temp * Xb_temp.t();
	    }

	    arma::mat D_sum (p_x*n_new, p_x, fill::zeros);
	    arma::mat D_temp = D[0];
	    D_sum(span(0, p_x - 1) , span::all) = D_temp;
	    for(int i = 2; i <= n_new ; i++){
	    	arma::mat D_temp = D[i-1];
	    	//Rcout <<  size(C_temp);
	    	D_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = D_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + D_temp; 
	    }
	    arma::mat D_sum_new (p_x*n_new, p_x, fill::zeros);
	    D_sum_new( span(  0, p_x - 1 ), span::all) =  D_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
	    for(int i = 2; i <= n_new ; i++){
	    	D_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = D_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - D_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
	    }

	    arma::mat D_sum_new_inv (p_x*n_new,p_x, fill::zeros);

	    for(int i = 1; i <= n_new ; i++){
	        vec eigval =  eig_sym(D_sum_new ( span((i-1)*p_x, (i*p_x)-1), span::all)  );
	        double add_pd =0.0;
	        double min_eigen = eigval(0);
	        if(min_eigen <= 0){
	            Rprintf("Invertiable! adding noise!");
	            if(min_eigen < 0){
	            	add_pd = (10)*fabs( min_eigen);
	            }
	            else
	            	add_pd = pow(10.0, -6.0);            
	        }
	        arma::mat noise(p_x,p_x, fill::eye);
	        D_sum_new_inv(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = inv( D_sum_new( span((i-1)*p_x, (i*p_x) -1  ), span::all  ) +  add_pd*noise );
	    }


	    
	    arma::mat phi_new (p_y, p_x*n_new, fill::zeros);


        arma::mat phi_compare = phi_hat;

        for(int i = 1; i <= n_new ; i++){
        	List E(n_new);
        	for(int j = 1; j <= n_new ; j++){       		
        		E[j-1] = D_sum_new( span(  ( std::max( j,i )-1)*p_x  , std::max(j,i )*p_x  -1 ), span::all)  * phi_hat(span::all , span( (j-1)*p_x, j*p_x -1  )).t() ;    			
        	}
        	

        	arma::mat E_new = E[0];
    		for ( int g = 1; g < n_new; g++ ) {
    			arma::mat E_temp = E[g];
    			E_new = E_new + E_temp;
    		}
    		

    		E_new =  E_new - D_sum_new(  span( (i-1)*p_x,  i*p_x -1  ), span::all)  * phi_hat( span::all, span( (i-1)*p_x, i*p_x -1) ).t();

    		arma::mat S =  C_sum_new(  span(  (i-1)*p_x , i*p_x -1 ), span::all )  - E_new;
   		
    		S = soft_full(S, lambda);

    		arma::mat phi_temp = D_sum_new_inv( span ((i-1)*p_x , i*p_x-1), span::all )  *  S;
    		
    		phi_temp = phi_temp.t();


    		phi_hat( span::all ,   span( (i-1)*p_x, i*p_x -1)  ) = phi_temp;
            phi_new( span::all ,   span( (i-1)*p_x, i*p_x -1)  ) = phi_temp;
                                                                  

        }

        // Rcout << phi_hat;

        arma::mat phi_temp_soft(p_y,p_x*n_new, fill::zeros);
        phi_temp_soft( span::all , span(0, p_x-1) ) = soft_full(phi_hat( span::all ,  span(0, p_x-1) ),lambda2);
        arma::mat temp_1 = phi_hat( span::all ,  span(0, p_x-1) );
        for(int z_1 = 2; z_1 <= n_new; z_1 ++){
            arma::mat temp_2 = temp_1 + phi_hat( span::all ,  span((z_1-1)*p_x, z_1*p_x-1) );
            arma::mat temp_1_soft = soft_full(temp_1,lambda2);
            arma::mat temp_2_soft = soft_full(temp_2,lambda2);
            phi_temp_soft(   span::all ,  span( (z_1-1)*p_x , z_1*p_x  -1 ) )= temp_2_soft - temp_1_soft;
            temp_1 = temp_2;
        }
        phi_new  = phi_temp_soft;


        // Rcout << phi_new;
        arma::mat data_y_m_2(T, p_y);
    	arma::mat beta_hat(p_y,p_x*n_new, fill::zeros);
    	beta_hat(span::all, span( 0, p_x -1)) =  phi_new(span::all, span( 0, p_x -1));
    	data_y_m_2(span(blocks[0]-1,  blocks[1]-2 ), span::all)=  data_y_m_all(span(blocks[0]-1,  blocks[1]-2 ), span::all) - data_x_m(span(blocks[0]-1,  blocks[1]-2 ), span::all)* beta_hat(span::all, span( 0, p_x -1)).t();

    	for(int i = 1; i <n_new; i++){
    		beta_hat(span::all, span( i*p_x, (i+1)*p_x -1)) =  phi_new(span::all, span( i*p_x, (i+1)*p_x -1)) + beta_hat(span::all, span( (i-1)*p_x, i*p_x -1));
    		data_y_m_2(span(blocks[i]-1,  blocks[i+1]-2 ), span::all)=  data_y_m_all(span(blocks[i]-1,  blocks[i+1]-2 ), span::all) - data_x_m(span(blocks[i]-1,  blocks[i+1]-2 ), span::all)* beta_hat(span::all, span( i*p_x, (i+1)*p_x -1)).t();
    	}


    	
    	arma::mat coef = solve(data_x_m_2, data_y_m_2); 
    	phi_hat_2  = coef.t();

        // mat abs_temp_2 = abs(phi_hat_2);
        // for(int i = 0; i < p_y; i++){
        //     for( int j = 0 ; j < fixed_l; j++){
        //         if ( abs_temp_2(i, j) > 0.05) {
        //             if( phi_hat_2(i, j ) > 0){
        //                 phi_hat_2(i, j ) = 0.05; 
        //             }
        //             if( phi_hat_2(i, j ) < 0){
        //                 phi_hat_2(i, j ) = -0.05; 
        //             }                    
        //         }
        //     }
        // }
        
        



        arma::mat abs_temp = abs(phi_new - phi_compare);
        double max_temp = abs_temp.max();
        if ( max_temp < tol) {
            phi_hat = phi_new; 
            break;
        } 
        if ( max_temp > tol ) {
        	phi_hat = phi_new; 
            // Rprintf( "%f \n", max_temp);

        }
        if (  max_temp > pow(10.0, 5.0)) {
            flag = 1;
            break;
        }

                
	}
	

    return List::create(Named("phi.hat")= phi_hat, Named("phi.hat.2")= phi_hat_2, Named("flag")= flag);	
}

// [[Rcpp::export]]
List var_break_fit_block(NumericMatrix data, double lambda, double lambda2, int q, int max_iteration, double tol , NumericMatrix initial_phi, NumericVector blocks, NumericVector cv_index ){

    int p = data.ncol(); int T = data.nrow(); int n_new = blocks.size() - 1;
    int n = T - q;


    arma::mat data_m(data.begin(), T, p);
    List Y_b(n_new);
    List y_b(n_new);
    //Rcout  << size(data_m)[1];

    for(int i = 0; i < n_new; i++) {
        y_b[i] = data_m( span(blocks[i]-1,  blocks[i+1]-2 ), span::all) ;
    }
    y_b[0] = data_m(  span(q, (blocks[1]-2) ), span::all );


    arma::mat Y( p*q, n);
    for(int i = (q-1); i < (T-1); i++) {
        for(int j = 1; j <= q; j++){
            Y(span((j-1)*p, j*p-1), i-(q-1) ) = data_m( i-j+1, span::all  ).t();
            //Y.submat( (j-1)*p,(i-q+1), j*p-1 ,(i-q+1) ) = data_m.submat(i-j+1, 0, i-j+1,  p-1 ).t();
        }       

    }
    
    Y_b[0] = Y( span(0 ,p*q-1) ,span(0, blocks[1]-q-2 ) );
    for(int i = 1; i < n_new; i++) {
        Y_b[i] = Y( span(0, p*q-1 ), span(blocks[i]-q-1 ,blocks[i+1]-q-2  )) ;
    }
    
    //NumericMatrix temp = y_b[0]; 
    //Rcout << temp.nrow();

    int cv_l = cv_index.size();
    if( cv_l >0){
        for(int t_1 = 0; t_1 < cv_l; t_1 ++){
            mat yb_temp = y_b[cv_index[t_1]-1]; 
            mat Yb_temp = Y_b[cv_index[t_1]-1]; 
            int tt =  yb_temp.n_rows;
            //Rcout << tt1;
            y_b[cv_index[t_1]-1] = yb_temp(span(0,tt-2), span::all );
            Y_b[cv_index[t_1]-1] = Yb_temp( span::all ,span(0,tt-2));

        }
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
        arma::mat yb_temp = y_b[j];
        arma::mat Yb_temp = Y_b[j]; 
        C[j] = Yb_temp * yb_temp;
    }  
    
    arma::mat C_sum (p*q*n_new, p, fill::zeros);
    //Rcout << C_sum(1,1);
    arma::mat C_temp = C[0];
    C_sum(span(0, p*q - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
        arma::mat C_temp = C[i-1];
        //Rcout <<  size(C_temp);
        C_sum(  span((i-1)*p*q, (i*p*q) -1  ), span::all ) = C_sum( span( (i-2)*p*q,(i-1)*p*q -1 ), span::all) + C_temp; 
    }
    arma::mat C_sum_new (p*q*n_new,p, fill::zeros);
    C_sum_new( span(   0, p*q - 1 ), span::all) =  C_sum(   span(  (n_new-1)*p*q,  n_new*p*q - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
        C_sum_new(  span((i-1)*p*q, (i*p*q) -1  ), span::all )  = C_sum(   span(  (n_new-1)*p*q,  n_new*p*q - 1  ), span::all ) - C_sum( span ((i-2)*p*q, (i-1)*p*q -1 ), span::all);
    }


    List D(n_new);
    for(int j = 0; j < n_new; j++){
        arma::mat Yb_temp = Y_b[j]; 
        D[j] = Yb_temp * Yb_temp.t();
    }

    arma::mat D_sum (p*q*n_new,p*q, fill::zeros);
    arma::mat D_temp = D[0];
    D_sum(span(0, p*q - 1) , span::all) = D_temp;
    for(int i = 2; i <= n_new ; i++){
        arma::mat D_temp = D[i-1];
        //Rcout <<  size(C_temp);
        D_sum(  span((i-1)*p*q, (i*p*q) -1  ), span::all ) = D_sum( span( (i-2)*p*q,(i-1)*p*q -1 ), span::all) + D_temp; 
    }
    arma::mat D_sum_new (p*q*n_new,p*q, fill::zeros);
    D_sum_new( span(  0, p*q - 1 ), span::all) =  D_sum(   span(  (n_new-1)*p*q,  n_new*p*q - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
        D_sum_new(  span((i-1)*p*q, (i*p*q) -1  ), span::all )  = D_sum(   span(  (n_new-1)*p*q,  n_new*p*q - 1  ), span::all ) - D_sum( span ((i-2)*p*q, (i-1)*p*q -1 ), span::all);
    }

    arma::mat D_sum_new_inv (p*q*n_new,p*q, fill::zeros);
    // for(int i = 1; i <= n_new ; i++){
    //  arma::mat noise(p*q,p*q, fill::eye);
    //  D_sum_new_inv(  span((i-1)*p*q, (i*p*q) -1  ), span::all )  = inv( D_sum_new( span((i-1)*p*q, (i*p*q) -1  ), span::all  ) +  pow(10, -6)*noise );
    // }

    for(int i = 1; i <= n_new ; i++){
        vec eigval =  eig_sym(D_sum_new ( span((i-1)*p*q, (i*p*q)-1), span::all)  );
        double add_pd =0.0;
        double min_eigen = eigval(0);
        //Rcout << eigval << '\n';
        //Rcout << min_eigen << '\n';
        if(min_eigen <= 0){
            Rprintf("Invertiable! adding noise!");
            if(min_eigen < 0){
            	add_pd = (10)*fabs( min_eigen);
            }
            else
            	add_pd = pow(10.0, -6.0); 
        }
        arma::mat noise(p*q,p*q, fill::eye);
        D_sum_new_inv(  span((i-1)*p*q, (i*p*q) -1  ), span::all )  = inv( D_sum_new( span((i-1)*p*q, (i*p*q) -1  ), span::all  ) +  add_pd*noise );
    }

    // mat D_new (n_new, q*p, fill::zeros);
    // for(int i = 0; i < n_new ; i++){
    //  mat Yb_temp = Y_b[i]; 
    //  mat D_new_temp  = diagvec(Yb_temp * Yb_temp.t() ); 
    //  D_new(i , span::all) =  reshape(D_new_temp, 1, q*p);
    // }


    arma::mat phi_hat(initial_phi.begin(), p, p*q*n_new);
    arma::mat phi_new (p, p*q*n_new, fill::zeros);

    int flag = 0;
    int l = 2;
    //Rcout <<  floor(0.5 * max_iteration); 
    //Rcout <<  (l == floor(0.5 * max_iteration)); 
    
    while( l < max_iteration){
        if(  l == floor(0.5 * max_iteration)  ){
            tol = (2)*tol;
        }
        if( l  == floor( 0.75 * max_iteration )) {
            tol = (4/2)*tol;
        }
        l = l+1; 
        arma::mat phi_compare = phi_hat;

        for(int i = 1; i <= n_new ; i++){
            List E(n_new);
            for(int j = 1; j <= n_new ; j++){
                
                E[j-1] = D_sum_new( span(  ( std::max( j,i )-1)*p*q  , std::max(j,i )*p*q  -1 ), span::all)  * phi_hat(span::all , span( (j-1)*p*q, j*p*q -1  )).t() ;              
            }
            

            arma::mat E_new = E[0];
            for ( int g = 1; g < n_new; g++ ) {
                arma::mat E_temp = E[g];
                E_new = E_new + E_temp;
            }


            

            E_new =  E_new - D_sum_new(  span( (i-1)*p*q,  i*p*q -1  ), span::all)  * phi_hat( span::all, span( (i-1)*p*q, i*p*q -1) ).t();

            arma::mat S =  C_sum_new(  span(  (i-1)*p*q , i*p*q -1 ), span::all )  - E_new;

            

            S = soft_full(S, lambda);

            arma::mat phi_temp = D_sum_new_inv( span ((i-1)*p*q , i*p*q-1), span::all )  *  S;
            
            phi_temp = phi_temp.t();


            phi_hat( span::all ,   span( (i-1)*p*q, i*p*q -1)  ) = phi_temp;
            phi_new( span::all ,   span( (i-1)*p*q, i*p*q -1)  ) = phi_temp;                                                                  

        }



        arma::mat phi_temp_soft(p,p*q*n_new, fill::zeros);
        phi_temp_soft( span::all , span(0, p*q-1) ) = soft_full(phi_hat( span::all ,  span(0, p*q-1) ),lambda2);
        arma::mat temp_1 = phi_hat( span::all ,  span(0, p*q-1) );
        for(int z_1 = 2; z_1 <= n_new; z_1 ++){
            arma::mat temp_2 = temp_1 + phi_hat( span::all ,  span((z_1-1)*p*q, z_1*p*q-1) );
            arma::mat temp_1_soft = soft_full(temp_1,lambda2);
            arma::mat temp_2_soft = soft_full(temp_2,lambda2);
            phi_temp_soft(   span::all ,  span( (z_1-1)*p*q , z_1*p*q  -1 ) )= temp_2_soft - temp_1_soft;
            temp_1 = temp_2;
        }
        phi_new  = phi_temp_soft;



        arma::mat abs_temp = abs(phi_new - phi_compare);
        //Rcout << max_temp;
        double max_temp = abs_temp.max();
        //Rcout << max_temp;
        //Rcout << '\n';
        if ( max_temp < tol) {
            break;
        } 
        if ( max_temp > tol ) {
            phi_hat = phi_new;    
            // Rprintf( "%f \n", max_temp);
        }
        if (  max_temp > pow(10.0, 5.0)) {
            Rprintf("NOT CONVERGED");
            flag = 1;
            break;
        }
                
    }
    

    return List::create(Named("phi.hat")= phi_hat,Named("flag")= flag);

    
}

//' lambda warm up for linear regression model
//'
//' @param data_y input matrix Y
//' @param data_x input matrix X 
//' @param blocks the vector of blocks 
//' @param cv_index the vector of indices for validation
//' @return a value for parameter lambda
// [[Rcpp::export]]
List lambda_warm_up_lm(NumericMatrix data_y, NumericMatrix data_x, NumericVector blocks, NumericVector cv_index ){

    
    int p_y = data_y.ncol(); int p_x = data_x.ncol(); int T = data_y.nrow(); 
    int n_new = blocks.size() - 1;
    int n = T;


    arma::mat data_y_m(data_y.begin(), T, p_y);
    arma::mat data_x_m(data_x.begin(), T, p_x);
    List y_b(n_new);

    for(int i = 0; i < n_new; i++) {
        y_b[i] = data_y_m( span(blocks[i]-1,  blocks[i+1]-2 ), span::all) ;
    }

    
    arma::mat X( p_x, n);
    for(int i = 0; i < T; i++) {
            X(span::all, i )= data_x_m( i , span::all ).t();
    }

    List X_b(n_new);   
    for(int i = 0; i < n_new; i++) {
        X_b[i] = X(span::all, span(blocks[i]-1, blocks[i+1]-2 ));
    }

    // Rcout <<  X(span::all, span(blocks[0]-1, blocks[1]-2 ));

    int cv_l = cv_index.size();
    if( cv_l > 0){
        for(int t_1 = 0; t_1 < cv_l; t_1 ++){
            arma::mat yb_temp = y_b[cv_index[t_1]-1]; 
            arma::mat Xb_temp = X_b[cv_index[t_1]-1]; 
            int tt =  yb_temp.n_rows;
            if( tt > 1 ){
                y_b[cv_index[t_1]-1] = yb_temp(span(0,tt-2), span::all );
                X_b[cv_index[t_1]-1] = Xb_temp( span::all ,span(0,tt-2));
            }
            else{
                Rcout << "Wrong!";
            }          

        }
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
        arma::mat yb_temp = y_b[j];
        arma::mat Xb_temp = X_b[j]; 
        C[j] = Xb_temp * yb_temp;
    }  
    
    arma::mat C_sum (p_x*n_new, p_y, fill::zeros);
    //Rcout << C_sum(1,1);
    arma::mat C_temp = C[0];
    C_sum(span(0, p_x - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
        arma::mat C_temp = C[i-1];
        //Rcout <<  size(C_temp);
        C_sum(  span((i-1)*p_x, (i*p_x) -1  ), span::all ) = C_sum( span( (i-2)*p_x,(i-1)*p_x -1 ), span::all) + C_temp; 
    }

    arma::mat C_sum_new (p_x*n_new, p_y, fill::zeros);
    C_sum_new( span(   0, p_x - 1 ), span::all) =  C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
        C_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )  = C_sum(   span(  (n_new-1)*p_x,  n_new*p_x - 1  ), span::all ) - C_sum( span ((i-2)*p_x, (i-1)*p_x -1 ), span::all);
    }

    double lambda_max = 0;
    for(int i = 1; i <= n_new; i++){
        double lambda_temp = abs(C_sum_new(  span((i-1)*p_x, (i*p_x) -1  ), span::all )).max(); 
        lambda_max = std::max(lambda_max,lambda_temp);

    }


    return List::create(Named("lambda_1_max")= lambda_max);
    
}

// [[Rcpp::export]]
List lambda_warm_up_var(NumericMatrix data, int q, NumericVector blocks, NumericVector cv_index ){

    int p = data.ncol(); int T = data.nrow(); int n_new = blocks.size() - 1;
    int n = T - q;


    arma::mat data_m(data.begin(), T, p);
    List Y_b(n_new);
    List y_b(n_new);
    //Rcout  << size(data_m)[1];

    for(int i = 0; i < n_new; i++) {
        y_b[i] = data_m( span(blocks[i]-1,  blocks[i+1]-2 ), span::all) ;
    }
    y_b[0] = data_m(  span(q, (blocks[1]-2) ), span::all );


    arma::mat Y( p*q, n);
    for(int i = (q-1); i < (T-1); i++) {
        for(int j = 1; j <= q; j++){
            Y(span((j-1)*p, j*p-1), i-(q-1) ) = data_m( i-j+1, span::all  ).t();
            //Y.submat( (j-1)*p,(i-q+1), j*p-1 ,(i-q+1) ) = data_m.submat(i-j+1, 0, i-j+1,  p-1 ).t();
        }       

    }
    
    Y_b[0] = Y( span(0 ,p*q-1) ,span(0, blocks[1]-q-2 ) );
    for(int i = 1; i < n_new; i++) {
        Y_b[i] = Y( span(0, p*q-1 ), span(blocks[i]-q-1 ,blocks[i+1]-q-2  )) ;
    }
    
    //NumericMatrix temp = y_b[0]; 
    //Rcout << temp.nrow();

    int cv_l = cv_index.size();
    if( cv_l >0){
        for(int t_1 = 0; t_1 < cv_l; t_1 ++){
            arma::mat yb_temp = y_b[cv_index[t_1]-1]; 
            arma::mat Yb_temp = Y_b[cv_index[t_1]-1]; 
            int tt =  yb_temp.n_rows;
            //Rcout << tt1;
            y_b[cv_index[t_1]-1] = yb_temp(span(0,tt-2), span::all );
            Y_b[cv_index[t_1]-1] = Yb_temp( span::all ,span(0,tt-2));

        }
    } 

    List C(n_new);
    for(int j = 0; j < n_new; j++){
        arma::mat yb_temp = y_b[j];
        arma::mat Yb_temp = Y_b[j]; 
        C[j] = Yb_temp * yb_temp;
    }  
    
    arma::mat C_sum (p*q*n_new, p, fill::zeros);
    //Rcout << C_sum(1,1);
    arma::mat C_temp = C[0];
    C_sum(span(0, p*q - 1) , span::all) = C_temp;
    for(int i = 2; i <= n_new ; i++){
        arma::mat C_temp = C[i-1];
        //Rcout <<  size(C_temp);
        C_sum(  span((i-1)*p*q, (i*p*q) -1  ), span::all ) = C_sum( span( (i-2)*p*q,(i-1)*p*q -1 ), span::all) + C_temp; 
    }
    arma::mat C_sum_new (p*q*n_new,p, fill::zeros);
    C_sum_new( span(   0, p*q - 1 ), span::all) =  C_sum(   span(  (n_new-1)*p*q,  n_new*p*q - 1  ), span::all );
    for(int i = 2; i <= n_new ; i++){
        C_sum_new(  span((i-1)*p*q, (i*p*q) -1  ), span::all )  = C_sum(   span(  (n_new-1)*p*q,  n_new*p*q - 1  ), span::all ) - C_sum( span ((i-2)*p*q, (i-1)*p*q -1 ), span::all);
    }


    double lambda_max = 0;
    for(int i = 1; i <= n_new; i++){
        double lambda_temp = abs(C_sum_new(  span((i-1)*p*q, (i*p*q) -1  ), span::all )).max(); 
        lambda_max = std::max(lambda_max,lambda_temp);

    }


    return List::create(Named("lambda_1_max")= lambda_max);

    
}

