#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fq_arith.h"
#include "monomial_mat.h"
#include "timing_and_stat.h"
#include "codes.h"
#include "LESS.h"
#include "rng.h"

#define GRN "\e[0;32m"
#define WHT "\e[0;37m"

#ifdef N_pad
#define NN N_pad
#else
#define NN N
#endif

#define CASE 1
#define NUM_TEST_ITERATIONS 100
#define PROB 0.05f

int sample_binom_p(){
    unsigned long M = 1<<20;
    unsigned long MAX = (unsigned long)(M/PROB);
    unsigned long p;
    randombytes((unsigned char*)&p, 8);
    p = p%MAX;
    if(p < M)
        return 1;
    else
        return 0;
}

/* Exhaustive testing of inverses mod Q */
void inverse_mod_tester(){
    uint32_t value[Q-1];
    uint32_t inverse[Q-1];
    for(uint32_t i=1; i <= Q-1; i++){
        value[i-1] = i;
        inverse[i-1] = fq_inv(i);
    }
    int all_ok = 1;
    for(uint32_t i=1; i <= Q-1; i++){
        if((value[i-1]*inverse[i-1]) % Q !=1){
           printf("%u*%u=%u\n",
                  value[i-1],
                  inverse[i-1],
                  (value[i-1]*inverse[i-1])%Q);
           all_ok = 0;
        }
    }
    if (all_ok){
        puts("All inverses on F_q ok");
    }
}


/* Test monomial matrix multiplication and inversion testing if
 * M*M^-1 = I, where I is a random monomial matrix */
void monomial_tester(){
    monomial_t mat1, mat2, id,idcheck;
    monomial_mat_id(&idcheck);
    monomial_mat_rnd(&mat1);
    monomial_mat_inv(&mat2,&mat1);
    monomial_mat_mul(&id,&mat2,&mat1);
    if( memcmp( &id,&idcheck,sizeof(monomial_t)) !=0 ){
       monomial_mat_pretty_print_name("M1",&mat1);
       monomial_mat_pretty_print_name("M1^-1",&mat2);
       monomial_mat_pretty_print_name("M1^-1 * M1",&id);
    } else {
        printf("Monomial arith test: ok\n");
    }
}



/*Generate a random full-rank G, keep multiplying G by a monomial and RREF'ing */
#define NUMBER_OF_GE_TESTS 10000
void gausselim_tester(){
    generator_mat_t G,GMul;
    uint8_t is_pivot_column[NN];

    /* randomly generate a non-singular G */
    do {
        generator_rnd(&G);
        memset(is_pivot_column,0,sizeof(is_pivot_column));
    } while ( generator_RREF(&G,is_pivot_column) == 0);

    /* Stress-test GE by repeatedly bringing in SF GQ*/
    monomial_t mat1;
    int full_rank,rref_ok, all_ok = 1;
    for(int i=0; (i<NUMBER_OF_GE_TESTS) && all_ok; i++ ){
        monomial_mat_rnd(&mat1);
        generator_monomial_mul(&GMul,&G,&mat1);
        memcpy(&G,&GMul,sizeof(generator_mat_t));
        memset(is_pivot_column,0,sizeof(is_pivot_column));
        full_rank = generator_RREF(&GMul,is_pivot_column);
        if(!full_rank){
            all_ok = 0;
            fprintf(stderr,"Singular Matrix (iter:%d)\n",i);
            generator_pretty_print_name("Pre-GE",&G);
            generator_pretty_print_name("Post-GE",&GMul);
        }

        /* check if the matrix is in reduced row echelon form i.e,
         *  - there are k pivots
         *  - each pivot appears as the first element of a row
         *  - is_pivot_column indicator array is correct
         */
        rref_ok = 1;
        for(int row_idx = 0; row_idx < K; row_idx++){
            int found_pivot_column = 0;
            while ( (GMul.values[row_idx][found_pivot_column] == 0) &&
                   (found_pivot_column < N) ){
                found_pivot_column++;
            }
            if ( (GMul.values[row_idx][found_pivot_column] != 1) ){
                fprintf(stderr,"row %d Pivot actually equal to %d\n",row_idx, GMul.values[row_idx][found_pivot_column]);
                     rref_ok = 0;
                }
            if ( (found_pivot_column >= N) ){
                fprintf(stderr,"row %d Pivot missing\n",row_idx);
                     rref_ok = 0;
                }
            if ( (is_pivot_column[found_pivot_column] != 1)){
                fprintf(stderr,"row %d indicator array mismatch\n",row_idx);
                     rref_ok = 0;
                }
        }

        if(full_rank && !rref_ok){
            fprintf(stderr,"RREF incorrect (iter:%d)\n",i);
            fprintf(stderr,"Pre-RREF\n");
            generator_pretty_print_name("Pre-RREF",&G);
            fprintf(stderr,"is_pivot = \n [ ");
            for(int x=0;x < N ;x++){fprintf(stderr," %d ",is_pivot_column[x]); }
            fprintf(stderr,"]\n");
            fprintf(stderr,"Post-RREF\n");
            generator_pretty_print_name("Post-RREF",&GMul);
            all_ok = 0;
        }

    }
    if(all_ok) {
        printf("GE test: ok\n");
    }
}


/* tests if G*M1*(M1^-1) = G*/
void gen_by_monom_tester(){
     generator_mat_t G = {0}, G2, Gcheck;
     uint8_t is_pivot_column[NN];
     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);
     monomial_t mat1, mat2;
     monomial_mat_rnd(&mat1);
     monomial_mat_inv(&mat2,&mat1);
     generator_monomial_mul(&G2,&G,&mat1);
     generator_monomial_mul(&Gcheck,&G2,&mat2);
     if( memcmp( &Gcheck,&G,sizeof(generator_mat_t)) !=0 ){
        generator_pretty_print_name("G",&G);
        generator_pretty_print_name("G*Q",&G2);
        generator_pretty_print_name("G*Q*Q^-1",&Gcheck);
     } else {
         printf("Generator-monomial multiplication: ok\n");
     }
}


/* draw random full rank G and pack it */
/* compute mu = Q_a^-1 Q_b */
/* compute G2 = G Q_a */
/* compute G3 = G Q_b */
/* compute Gcheck = G2 mu */

void rref_gen_by_monom_tester(){
     generator_mat_t G, G2, G3, Gcheck;

     monomial_t Q_a,Q_a_inv,Q_b,mu,Qcheck;
     monomial_mat_rnd(&Q_a);
     monomial_mat_inv(&Q_a_inv,&Q_a);
     monomial_mat_rnd(&Q_b);

     uint8_t is_pivot_column[NN];
     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);

     generator_monomial_mul(&G2,&G,&Q_a);
     if (generator_RREF(&G2,is_pivot_column) != 1){
         printf("G2=G Q_a: singular\n");

     };
     generator_monomial_mul(&G3,&G,&Q_b);
     if (generator_RREF(&G3,is_pivot_column) != 1){
         printf("G3=G Q_b: singular\n");

     };
     monomial_mat_mul(&mu,&Q_a_inv,&Q_b);


    monomial_mat_mul(&Qcheck,&Q_a,&mu);
    if( memcmp( &Q_b,&Qcheck,sizeof(monomial_t)) !=0 ){
        monomial_mat_pretty_print_name("mu",&mu);
        monomial_mat_pretty_print_name("Q_a",&Q_a);
        monomial_mat_pretty_print_name("Qcheck",&Qcheck);
        monomial_mat_pretty_print_name("Q_b",&Q_b);
        fprintf(stderr,"Q_a mu != Q_b\n");
    }


     generator_monomial_mul(&Gcheck,&G2,&mu);
     generator_RREF(&Gcheck,is_pivot_column);
     if (generator_RREF(&Gcheck,is_pivot_column) != 1){
         printf("Gcheck=G2 mu: singular\n");

     };

     if( memcmp( &Gcheck,&G3,sizeof(generator_mat_t)) !=0 ){
         printf("SF Generator-monomial multiplication: not ok\n");
        generator_pretty_print_name("G",&G);
        generator_pretty_print_name("G2",&G2);
        generator_pretty_print_name("G3",&G3);
        generator_pretty_print_name("Gcheck",&Gcheck);
        monomial_mat_print_exp_name("Q_a",&Q_a);
        monomial_mat_print_exp_name("Q_a_inv",&Q_a_inv);
        monomial_mat_print_exp_name("Q_b",&Q_b);
        monomial_mat_print_exp_name("mu",&mu);
        monomial_mat_print_exp_name("Qcheck",&Qcheck);
     } else {
         printf("SF Generator-monomial multiplication: ok\n");
     }
}

void rref_gen_compress_tester(){
     generator_mat_t G = {0}, Gcheck;
     rref_generator_mat_t SF_G;
     uint8_t is_pivot_column[NN];

     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);

     memcpy(&Gcheck,&G, sizeof(G));
     generator_rref_compact(&SF_G,&G,is_pivot_column);
     generator_rnd(&G); /* fill with garbage to elicit faults */
     generator_rref_expand(&G,&SF_G);

    if( memcmp( &Gcheck,&G,sizeof(generator_mat_t)) !=0 ){
        printf("Generator SF compression: ko\n");
       fprintf(stderr," Comp-decomp\n");
       generator_pretty_print_name("G",&G);
       fprintf(stderr,"is_pivot = \n [ ");
       for(int x=0;x < N ;x++){fprintf(stderr," %d ",is_pivot_column[x]); }
       fprintf(stderr,"]\n");

       generator_rref_pretty_print_name("RREF-G",&SF_G);

       fprintf(stderr," Reference\n");
       generator_pretty_print_name("Gcheck",&Gcheck);
    } else {
        printf("Generator SF compression: ok\n");
    }
}

void mono_is_compress_tester(){
    monomial_action_IS_t Q_a, Qcheck;
    uint8_t compressed [MONO_ACTION_PACKEDBYTES];

    monomial_t mono_rnd;
    monomial_mat_rnd(&mono_rnd);

    // Create random q
    for (int i = 0; i < K; i++) {
        Q_a.coefficients[i] = mono_rnd.coefficients[i];
        Q_a.permutation[i] = mono_rnd.permutation[i];
    }

     compress_monom_action(compressed,&Q_a);
     expand_to_monom_action(&Qcheck,compressed);

    if( memcmp( &Qcheck,&Q_a,sizeof(monomial_action_IS_t)) !=0 ){
        printf("Monomial Action compression: ko\n");

       fprintf(stderr,"perm = [");
       for(int i = 0; i < K-1; i++) {
          fprintf(stderr,"%03u, ",Q_a.permutation[i]);
       }
       fprintf(stderr,"%03u ]\n",Q_a.permutation[K-1]);
       fprintf(stderr,"coeffs = [");
       for(int i = 0; i < K-1; i++) {
          fprintf(stderr,"%03u, ",Q_a.coefficients[i]);
       }
       fprintf(stderr,"%03u ]\n",Q_a.coefficients[K-1]);

       fprintf(stderr,"\n\n\n");
       fprintf(stderr,"perm = [");
       for(int i = 0; i < K-1; i++
        ) {
          fprintf(stderr,"%03u, ",Qcheck.permutation[i]);
       }
       fprintf(stderr,"%03u ]\n",Qcheck.permutation[K-1]);
       fprintf(stderr,"coeffs = [");
       for(int i = 0; i < K-1; i++) {
          fprintf(stderr,"%03u, ",Qcheck.coefficients[i]);
       }
       fprintf(stderr,"%03u ]\n",Qcheck.coefficients[K-1]);

    } else {
        printf("Monomial Action compression: ok\n");
    }

}


void rref_gen_byte_compress_tester(){
     generator_mat_t G = {0}, Gcheck;
     uint8_t G_compressed [RREF_MAT_PACKEDBYTES];
     uint8_t is_pivot_column[NN];

     /* randomly generate a non-singular G */
     do {
         generator_rnd(&G);
         memset(is_pivot_column,0,sizeof(is_pivot_column));
     } while ( generator_RREF(&G,is_pivot_column) == 0);

     memcpy(&Gcheck,&G, sizeof(G));
     compress_rref(G_compressed,&G,is_pivot_column);
     generator_rnd(&G); /* fill with garbage to elicit faults */
     expand_to_rref(&G,G_compressed);

    if( memcmp( &Gcheck,&G,sizeof(generator_mat_t)) !=0 ){
        printf("Generator SF byte compression: ko\n");
       fprintf(stderr," Comp-decomp\n");
       generator_pretty_print_name("G",&G);

       fprintf(stderr,"is_pivot = \n [ ");
       for(int x=0;x < N ;x++){fprintf(stderr," %d ",is_pivot_column[x]); }
       fprintf(stderr,"]\n");

       fprintf(stderr," \n\n\n\n\n\n\n\n\nReference\n");
       generator_pretty_print_name("Gcheck",&Gcheck);
    } else {
        printf("Generator SF compression: ok\n");
    }
}

void info(){
    fprintf(stderr,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stderr,"num. keypairs = %d\n",NUM_KEYPAIRS);
    fprintf(stderr,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stderr,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stderr,"Public key %luB\n", sizeof(pubkey_t));
    fprintf(stderr,"Signature: %luB\n", sizeof(sig_t));

}

/* returns 1 if the test is successful, 0 otherwise */
int LESS_sign_verify_test(){
    pubkey_t pk;
    prikey_t sk;
    sig_t signature;
    char message[8] = "Signme!";
    LESS_keygen(&sk,&pk);
    LESS_sign(&sk,message,8,&signature);
    int is_signature_ok;
    is_signature_ok = LESS_verify(&pk,message,8,&signature);
    // fprintf(stderr,"Keygen-Sign-Verify: %s", is_signature_ok == 1 ? "functional\n": "not functional\n" );
    return is_signature_ok;
}


// *** new additional functions ****//

#define TO_PUBLISH 0
#define NOT_TO_PUBLISH 1
#define LEFT_CHILD(i) (2*i+1)
#define RIGHT_CHILD(i) (2*i+2)
#define PARENT(i) ((i-1)/2)
#define SIBLING(i) ( ((i)%2) ? i+1 : i-1 )
#define IS_LEFT_SIBLING(i) (i%2)

void find_partial_monomial_update2(monomial_t *Q_partial, monomial_action_IS_t *Q_in, monomial_action_IS_t *Q_out, int* IS_J){
    
    monomial_t Q_in_inv;
    
    for(int i=0; i<N; i++){
        Q_in_inv.permutation[i]=N+1;
        Q_in_inv.coefficients[i]=0;
    }
    for(int i=0; i<K; i++){
        Q_in_inv.permutation[Q_in->permutation[i]] = i;
        Q_in_inv.coefficients[Q_in->permutation[i]] = Q_in->coefficients[i];
    }
    
    for(int i=0;i<K;i++)
        IS_J[Q_in->permutation[i]] = 1;


    for(int i=0; i<N; i++){
        if(Q_in_inv.permutation[i]!=(N+1) && Q_partial->permutation[i]==(N+1)){
            
            Q_partial->permutation[i] = Q_out->permutation[Q_in_inv.permutation[i]];
            Q_partial->coefficients[i] = fq_red((FQ_DOUBLEPREC) Q_out->coefficients[Q_in_inv.permutation[i]] * fq_inv(Q_in_inv.coefficients[i]));
        }
    }

}/*Our function*/

void find_partial_monomial_update(monomial_t *Q_partial, monomial_action_IS_t *Q_in, monomial_action_IS_t *Q_out){
    
    monomial_t Q_in_inv;
    
    for(int i=0; i<N; i++){
        Q_in_inv.permutation[i]=N+1;
        Q_in_inv.coefficients[i]=0;
    }
    for(int i=0; i<K; i++){
        Q_in_inv.permutation[Q_in->permutation[i]] = i;
        Q_in_inv.coefficients[Q_in->permutation[i]] = Q_in->coefficients[i];
    }
    
    for(int i=0; i<N; i++){
        if(Q_in_inv.permutation[i]!=(N+1) && Q_partial->permutation[i]==(N+1)){
            
            Q_partial->permutation[i] = Q_out->permutation[Q_in_inv.permutation[i]];
            Q_partial->coefficients[i] = fq_red((FQ_DOUBLEPREC) Q_out->coefficients[Q_in_inv.permutation[i]] * fq_inv(Q_in_inv.coefficients[i]));
        }
    }

}/*Our function*/

void partial_monomial_initialize(monomial_t *Q_partial){
   for(int i=0; i<N; i++){
        Q_partial->permutation[i]=N+1;
        Q_partial->coefficients[i]=0;
    }
}/*Our function*/

void partial_monomial_transpose(monomial_t *Q_partial_transpose, monomial_t *Q_partial){
   for(int i=0; i<N; i++){
         Q_partial_transpose->permutation[i] = N+1;
         Q_partial_transpose->coefficients[i]= 0; 
      
    }
   for(int i=0; i<N; i++){
      if(Q_partial->permutation[i]!=N+1){
         Q_partial_transpose->permutation[Q_partial->permutation[i]] = i;
         Q_partial_transpose->coefficients[Q_partial->permutation[i]]= Q_partial->coefficients[i]; 
      }
    }
}/*Our function*/

int update_count_mono(monomial_t *Q_partial){
   int count = 0;
   for(int i=0; i<N; i++){
      if(Q_partial->permutation[i]!=(N+1))
         count++;      
    }
   return count;
}/*Our function*/

int Error_count(monomial_t *Q_partial, monomial_t *Q_transpose){
   int count = 0;
   for(int i=0; i<N; i++){
      if(Q_partial->permutation[i]!=(N+1))
         if(Q_partial->permutation[i]!=Q_transpose->permutation[i])
            return 1;      
    }
   return count;
}


int generate_partial_seed_tree(unsigned char
                                  seed_tree[NUM_NODES_OF_SEED_TREE *
                                                               SEED_LENGTH_BYTES],
                                  const unsigned char root_seed[SEED_LENGTH_BYTES],
                                  const unsigned char salt[HASH_DIGEST_LENGTH])
{
   /* input buffer to the CSPRNG, contains a salt, the seed to be expanded
    * and the integer index of the node being expanded for domain separation */
   const uint32_t csprng_input_len = HASH_DIGEST_LENGTH +
                                     SEED_LENGTH_BYTES +
                                     sizeof(uint32_t);
   unsigned char csprng_input[csprng_input_len];
   SHAKE_STATE_STRUCT tree_csprng_state;

   memcpy(csprng_input, salt, HASH_DIGEST_LENGTH);

   /* Set the root seed in the tree from the received parameter */
   uint32_t k=0,count=0;
   memcpy(seed_tree+SEED_LENGTH_BYTES,root_seed,SEED_LENGTH_BYTES);
   for (uint32_t i = 1; i < NUM_LEAVES_OF_SEED_TREE-1; ) {
      /* prepare the CSPRNG input to expand the children of node i */
      memcpy(csprng_input + HASH_DIGEST_LENGTH,
             seed_tree + i*SEED_LENGTH_BYTES,
             SEED_LENGTH_BYTES);
      *((uint32_t *)(csprng_input + HASH_DIGEST_LENGTH + SEED_LENGTH_BYTES)) = i;
      /* expand the children (stored contiguously) */
      initialize_csprng(&tree_csprng_state, csprng_input, csprng_input_len);
      csprng_randombytes(seed_tree + LEFT_CHILD(i)*SEED_LENGTH_BYTES,
                         2*SEED_LENGTH_BYTES,
                         &tree_csprng_state);
      count++;
      if(count==(1<<k)){
         k++;
         count=0;
         i=(uint32_t)((1<<(k+1))-1);
      }
      else{
         i++;
      }
   }
   return (1<<k);
} /*Generate partial ephem_seeds from the node tree[1]*//*Our function*/

/*************************************************
//                New functions
//                Updated versions
//
************************************************/
int compare(generator_mat_t *G, generator_mat_t *G_tilde, FQ_ELEM multiplier, int col){
    
    int i, j;
    for(i=0; i<N; i++){
        for(j=0; j<K; j++){
            FQ_ELEM tmp = fq_red((FQ_DOUBLEPREC)G_tilde->values[j][i] * (FQ_DOUBLEPREC)multiplier);
            if(tmp!=G->values[j][col])
                break;
        }
        if(j==K)
            return i;
    }
    return N+1;
}
void generator_monomial_mul_update(invertible_mat_t *res,
                            const generator_mat_t *const G,
                            const monomial_action_IS_t *const monom)
{
   for(int src_col_idx = 0; src_col_idx < K; src_col_idx++) {
      for(int row_idx = 0; row_idx < K; row_idx++) {
         res->values[row_idx][monom->permutation[src_col_idx]] =
            fq_red( (FQ_DOUBLEPREC) G->values[row_idx][src_col_idx] *
                    (FQ_DOUBLEPREC) monom->coefficients[src_col_idx] );
      }
   }
}

static inline
void swap_rows(FQ_ELEM r[N],FQ_ELEM s[N])
{
   FQ_ELEM tmp;
   for(int i=0; i<N; i++) {
      tmp = r[i];
      r[i] = s[i];
      s[i] = tmp;
   }
} /* end swap_rows */

int matrix_inverse(invertible_mat_t *G)
{
    generator_mat_t G_tilde;
    for(int row_to_reduce = 0; row_to_reduce < K; row_to_reduce++){
        for(int col_to_reduce = 0; col_to_reduce < K; col_to_reduce++){
            G_tilde.values[row_to_reduce][col_to_reduce] = G->values[row_to_reduce][col_to_reduce];
        }
        for(int col_to_reduce = K; col_to_reduce < N; col_to_reduce++){
            G_tilde.values[row_to_reduce][col_to_reduce] = 0;
            if(col_to_reduce== row_to_reduce+K)
                G_tilde.values[row_to_reduce][col_to_reduce] = 1;
        }
    }

   for(int row_to_reduce = 0; row_to_reduce < K; row_to_reduce++) {
      int pivot_row = row_to_reduce;
      /*start by searching the pivot in the col = row*/
      int pivot_column = row_to_reduce;
      while( (pivot_column < N) &&
             (G_tilde.values[pivot_row][pivot_column] == 0) ) {

         while ( (pivot_row < K) &&
                 (G_tilde.values[pivot_row][pivot_column] == 0) ) {
            pivot_row++;
         }
         if(pivot_row >= K) { /*entire column tail swept*/
            pivot_column++; /* move to next col */
            pivot_row = row_to_reduce; /*starting from row to red */
         }
      }
      if ( pivot_column >=N ) {
         return 0; /* no pivot candidates left, report failure */
      }
      //is_pivot_column[pivot_column] = 1; /* pivot found, mark the column*/


      /* if we found the pivot on a row which has an index > pivot_column
       * we need to swap the rows */
      if (row_to_reduce != pivot_row) {
         swap_rows(G_tilde.values[row_to_reduce],G_tilde.values[pivot_row]);
      }
      pivot_row = row_to_reduce; /* row with pivot now in place */


      /* Compute rescaling factor */
      FQ_DOUBLEPREC scaling_factor = fq_inv(G_tilde.values[pivot_row][pivot_column]);

      /* rescale pivot row to have pivot = 1. Values at the left of the pivot
       * are already set to zero by previous iterations */
      for(int i = pivot_column; i < N; i++) {
         G_tilde.values[pivot_row][i] = fq_red( (FQ_DOUBLEPREC) scaling_factor *
                                           (FQ_DOUBLEPREC) (G_tilde.values[pivot_row][i]) );
      }

      /* Subtract the now placed and reduced pivot rows, from the others,
       * after rescaling it */
      for(int row_idx = 0; row_idx < K; row_idx++) {
         if (row_idx != pivot_row) {
            FQ_DOUBLEPREC multiplier = G_tilde.values[row_idx][pivot_column];
            /* all elements before the pivot in the pivot row are null, no need to
             * subtract them from other rows. */
            for(int col_idx = 0; col_idx < N; col_idx++) {
               FQ_DOUBLEPREC tmp;
               tmp = fq_red( (FQ_DOUBLEPREC) multiplier *
                             (FQ_DOUBLEPREC) G_tilde.values[pivot_row][col_idx] );

               tmp = (FQ_DOUBLEPREC) Q + (FQ_DOUBLEPREC) G_tilde.values[row_idx][col_idx] - tmp;
               tmp = fq_red(tmp);

               G_tilde.values[row_idx][col_idx] = tmp;
            }
         }
      }
   }
   for(int row_to_reduce = 0; row_to_reduce < K; row_to_reduce++){
        for(int col_to_reduce = 0; col_to_reduce < K; col_to_reduce++){
            G->values[row_to_reduce][col_to_reduce] = G_tilde.values[row_to_reduce][K+col_to_reduce];

            if((row_to_reduce == col_to_reduce) && (G_tilde.values[row_to_reduce][col_to_reduce] != 1)) printf("\nError\n");
            if((row_to_reduce != col_to_reduce) && (G_tilde.values[row_to_reduce][col_to_reduce] != 0)) printf("\nError\n");
        }
    }


   return 1;
} /* end generator_RREF */

void mul_invertible_generator(generator_mat_t *G_tilde, invertible_mat_t *G, generator_mat_t *result_G){
    for(int i=0; i<K; i++){
        for(int j=0; j<N; j++){
            FQ_ELEM sum = 0;
            for(int r=0; r<K; r++){
                sum=fq_red((FQ_DOUBLEPREC)sum+fq_red((FQ_DOUBLEPREC)G->values[i][r]*(FQ_DOUBLEPREC)result_G->values[r][j]));
            }
            G_tilde->values[i][j] = sum;
        }
    }
}

void mul_invertible_invertible(invertible_mat_t *res, invertible_mat_t *G_1, invertible_mat_t *G_2){
    
    for(int i=0; i<K; i++){
        for(int j=0; j<K; j++){
            FQ_ELEM sum = 0;
            for(int r=0; r<K; r++){
                sum=fq_red((FQ_DOUBLEPREC)sum + fq_red((FQ_DOUBLEPREC)G_1->values[i][r]*(FQ_DOUBLEPREC)G_2->values[r][j]));
            }
            res->values[i][j] = sum;
        }
    }
}

int inverse(generator_mat_t G, monomial_action_IS_t Q_in, invertible_mat_t S){
    invertible_mat_t result_G;
    generator_monomial_mul_update(&result_G,
                             &G,
                             &Q_in);
    generator_inverse(&result_G);
    
}
int SIM_secret_key_single_fault(int *num_cols);
int SIM_recover_full_secret_key(int *num_sign);
int SIM_recover_full_secret_key_prob(int *num_sign);


/// End of new functions

#define failure 0
#define success 1


int main(int argc, char* argv[]){
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"012345678912345",
                      16);
    fprintf(stderr,"LESS reference implementation functional testbench\n");
    info();
    
    int count = 0;

    int check;

    //key_recovery_full(&count);  
    
    for(int i=0; i<NUM_TEST_ITERATIONS; i++){
        #if CASE==1
        
        check=SIM_secret_key_single_fault_new(&count);

        #elif CASE==2
        
        check=SIM_recover_full_secret_key_new(&count);

        #elif CASE==3

        check=SIM_recover_full_secret_key_prob_new(&count);

        #endif

        if(check == failure)
            printf("Failure at %d\n",i);    
    }
    #if CASE==1

    printf("\nThe average number of secret columns recovered with single fault: %f\n", (count)/(float)NUM_TEST_ITERATIONS);

    #elif CASE==2

    printf("\nThe average number of signatures required to recover full secret: %f\n", (count)/(float)NUM_TEST_ITERATIONS);

    #elif CASE==3

    printf("\nThe average number of signatures required to recover full secret: %f\n", (count)/(float)NUM_TEST_ITERATIONS);
    printf("Here the assumption is that the fault can be induced with probability %lf\n", PROB);

    #endif
    
    return 0;
}

/*******************************Updated Attack simulation function*****************/

int test_key_recovery(int* num_cols){  
    
    pubkey_t pk;
    prikey_t sk;
    LESS_keygen(&sk,&pk);

    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
    unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
    for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }
    monomial_t final_mono_secret[NUM_KEYPAIRS-1];
    monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];

    for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
    }

    monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

    for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
    }

    int count_mono_coeff[NUM_KEYPAIRS-1]={0};
    int total_count_coeff = 0;

    sig_t signature;
    char message[8] = "Signme!";      
    LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

    unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
    memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
    unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

    int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
    unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

    monomial_t Q_tilde;
    normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
    #ifdef COMPRESS_CMT_COLUMNS
    uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
    #endif
    monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
    generator_mat_t full_G0;
    generator_rref_expand(&full_G0,&G0_rref);

    int IS_J[N]={0}, IS_J_star={0};
    int employed_monoms = 0, F;
    monomial_action_IS_t mono_action;
    for(int i = 0; i < t; i++) {
        F = fixed_weight_string[i];
        if ( F!= 0){
           // if(count_mono_coeff[F-1]<N) {
               monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
               prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

               expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
               find_partial_monomial_update2(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action, IS_J);



               total_count_coeff-=count_mono_coeff[F-1];
               count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
               total_count_coeff+=count_mono_coeff[F-1];


               
            //}  
        employed_monoms++;
        break;
        }

    }

    monomial_t mono_partial_secret_inv;
    for(int i=0;i<N;i++)
    {
        mono_partial_secret_inv.permutation[i] = N+1;
        mono_partial_secret_inv.coefficients[i] = 0;
    }

    for(int i=0;i<N;i++)
    {
        if(mono_partial_secret[F-1].permutation[i] != N+1)
        {
            mono_partial_secret_inv.permutation[mono_partial_secret[F-1].permutation[i]] = i;
            mono_partial_secret_inv.coefficients[mono_partial_secret[F-1].permutation[i]] = mono_partial_secret[F-1].coefficients[i];
        }
    }

    for(int i=0;i<N;i++) 
                        printf("%d  ", mono_partial_secret[F-1].permutation[i]);
    printf("\n--------------------------\n");

    for(int i=0;i<N;i++) 
                        printf("%d  ", mono_partial_secret_inv.permutation[i]);
    printf("\n--------------------------\n");


    for(int i=0;i<N;i++) 
                        printf("%d  ", private_Q_inv[F-1].permutation[i]);
    printf("\n--------------------------\n");


    generator_mat_t G_hat;
    expand_to_rref(&G_hat, pk.SF_G[F-1]);

    invertible_mat_t G_hat_J_star, G_0_J, S, G_tmp;

    int k=0;
    /*
    for(int i=0; i<N;i++){
        if(IS_J[i] == 1){
            for(int j=0;j<K;j++){
                G_0_J.values[k][j] = full_G0.values[i][j];
            }
        }
        k++;
    }
    */
    for(int i=0; i<N;i++){
        if(mono_partial_secret[F-1].permutation[i] != N+1){
            for(int j=0;j<K;j++){
                FQ_ELEM tmp = fq_inv(mono_partial_secret[F-1].coefficients[i]);
                tmp = fq_red((FQ_DOUBLEPREC)tmp * (FQ_DOUBLEPREC)full_G0.values[j][i]);
                G_0_J.values[j][k] = tmp;
            }
            k++;
        }
    }

    k=0;
    for(int i=0;i<N;i++){
        if(mono_partial_secret[F-1].permutation[i] != N+1){
            for(int j=0;j<K;j++){
                G_hat_J_star.values[j][k] = G_hat.values[j][mono_partial_secret[F-1].permutation[i]];
            }
            k++;
        }
        
    }

    generator_mat_t G_left, G_right;

    G_tmp = G_0_J;
    matrix_inverse(&G_0_J);
    
    mul_invertible_invertible(&S, &G_hat_J_star, &G_0_J);
    
    /*
    monomial_mat_inv(&mono_partial_secret_inv, &private_Q_inv[F-1]);

    generator_monomial_mul(&G_left, &full_G0, &mono_partial_secret_inv);

    mul_invertible_generator(&G_right, &S, &G_left);

    for(int i=0;i<K;i++){
        for(int j=0;j<N;j++)
        {
            //if(j%12==0) printf("\n");
            if(G_right.values[i][j] != G_hat.values[i][j]) printf("1 ");
        }
        printf("\n------ ----------------- --------------------\n");
    }
    */
    
    matrix_inverse(&S);

    
    //mul_invertible_invertible(&G_0_J, &S, &G_hat_J_star);

    //for(int i=0;i<K;i++)
    //{
    //    for(int j=0;j<K;j++)
    //    {
    //        if(G_tmp.values[i][j] != G_0_J.values[i][j])
    //            printf("\nError\n");
    //    }
    //}
    mul_invertible_generator(&G_left, &S, &G_hat);

    //monomial_mat_inv(&mono_partial_secret_inv, &private_Q_inv[F-1]);

    //generator_monomial_mul(&G_right, &full_G0, &mono_partial_secret_inv);

    //for(int i=0;i<K;i++){
    //    for(int j=0;j<N;j++)
    //    {
            //if(j%12==0) printf("\n");
    //        if(G_right.values[i][j] != G_left.values[i][j]) printf("1 ");
    //    }
        //printf("\n------ ----------------- --------------------\n");
    //}
    

    for(int i=0;i<N;i++)
    {
        for(FQ_ELEM e = 1; e<Q; e++){
            int j = compare(&G_left, &full_G0, e, i);
            if(j!=N+1) printf("%d ", j);
        }
    }
    printf("\n");
    
    


    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        partial_monomial_transpose(&final_mono_secret[i], &private_Q_inv[i]);
    }

    int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(mono_partial_secret[i].permutation[j]!=N+1){
                if(mono_partial_secret[i].coefficients[j]!=final_mono_secret[i].coefficients[j]
                    ||mono_partial_secret[i].permutation[j]!=final_mono_secret[i].permutation[j]){
                        fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                    res = failure;
             }
            }
        }
    }
    res = success;

    if(res == success)
        *num_cols += total_count_coeff;
    return res;

}


int key_recovery_full(int* num_cols){  
    
    pubkey_t pk;
    prikey_t sk;
    LESS_keygen(&sk,&pk);

    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
    unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
    for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }
    monomial_t final_mono_secret[NUM_KEYPAIRS-1];
    monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];

    for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
    }

    monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

    for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
    }

    int count_mono_coeff[NUM_KEYPAIRS-1]={0};
    int total_count_coeff = 0;

    sig_t signature;
    char message[8] = "Signme!";      
    LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

    unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
    memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
    unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

    int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
    unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

    monomial_t Q_tilde;
    normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
    #ifdef COMPRESS_CMT_COLUMNS
    uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
    #endif
    monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
    generator_mat_t full_G0;
    generator_rref_expand(&full_G0,&G0_rref);

    int IS_J[N]={0}, IS_J_star={0};
    int employed_monoms = 0, F;
    monomial_action_IS_t mono_action;
    for(int i = 0; i < t; i++) {
        F = fixed_weight_string[i];
        if ( F!= 0){
            printf("%d\n",F);
            if(count_mono_coeff[F-1]<N) {
                monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
                prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

                expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
                find_partial_monomial_update2(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action, IS_J);



                total_count_coeff-=count_mono_coeff[F-1];
                
                total_count_coeff+=count_mono_coeff[F-1];

                generator_mat_t G_hat;
                expand_to_rref(&G_hat, pk.SF_G[F-1]);

                invertible_mat_t G_hat_J_star, G_0_J, S, G_tmp;

                int k=0;
                for(int x=0; x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            FQ_ELEM tmp = fq_inv(mono_partial_secret[F-1].coefficients[x]);
                            tmp = fq_red((FQ_DOUBLEPREC)tmp * (FQ_DOUBLEPREC)full_G0.values[j][x]);
                            G_0_J.values[j][k] = tmp;
                        }
                        k++;
                    }
                }

                k=0;
                for(int x=0;x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            G_hat_J_star.values[j][k] = G_hat.values[j][mono_partial_secret[F-1].permutation[x]];
                        }
                        k++;
                    }
                    
                }

                generator_mat_t G_left, G_right;

                G_tmp = G_0_J;

                matrix_inverse(&G_0_J);
                mul_invertible_invertible(&S, &G_hat_J_star, &G_0_J);
                matrix_inverse(&S);
                mul_invertible_generator(&G_left, &S, &G_hat);
                for(int x=0;x<N;x++)
                {
                    for(FQ_ELEM e = 1; e<Q; e++){
                        int j = compare(&G_left, &full_G0, e, x);
                        //if(j!=N+1) printf("%d %d\n", e, j);
                        if(j!=N+1) {
                            mono_partial_secret[F-1].permutation[x] = j;
                            mono_partial_secret[F-1].coefficients[x] = fq_inv(e);
                        }
                    }
                }
                count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
            } 

        employed_monoms++;
        
        }

    }

    //fprintf(stderr, "%d \n", fq_red((FQ_DOUBLEPREC)103*(FQ_DOUBLEPREC)37));

    //for(int i=0; i<NUM_KEYPAIRS-1; i++){
    //    partial_monomial_transpose(&final_mono_secret[i], &private_Q_inv[i]);
    //}

    printf("\n------------------\n");
    int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        res = success;
        for(int j=0; j<N; j++){
            
                if(mono_partial_secret[i].coefficients[j]!=private_Q_inv[i].coefficients[j]
                    ||mono_partial_secret[i].permutation[j]!=private_Q_inv[i].permutation[j]){
                       // fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                    res = failure;
             }
        }
        if(res == failure)
            printf("%d\n", i+1);
    }
    res = success;

    if(res == success)
        *num_cols += total_count_coeff;
    return res;

}

// Attack simulations //


int SIM_secret_key_single_fault(int *num_cols){  
    
    pubkey_t pk;
    prikey_t sk;
    LESS_keygen(&sk,&pk);

    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
    unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
    for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }
    monomial_t final_mono_secret[NUM_KEYPAIRS-1];
    monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];

    for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
    }

    monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

    for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
    }

    int count_mono_coeff[NUM_KEYPAIRS-1]={0};
    int total_count_coeff = 0;

    sig_t signature;
    char message[8] = "Signme!";      
    LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

    unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
    memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
    unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

    int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
    unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

    monomial_t Q_tilde;
    normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
    #ifdef COMPRESS_CMT_COLUMNS
    uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
    #endif
    monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
    generator_mat_t full_G0;
    generator_rref_expand(&full_G0,&G0_rref);

    int employed_monoms = 0;
    monomial_action_IS_t mono_action;
    for(int i = 0; i < t; i++) {
        int F = fixed_weight_string[i];
        if ( F!= 0){
            if(count_mono_coeff[F-1]<N) {
               monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
               prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

               expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
               find_partial_monomial_update(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action);
               total_count_coeff-=count_mono_coeff[F-1];
               count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
               total_count_coeff+=count_mono_coeff[F-1];
            }  
        employed_monoms++;
        }
    }


    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        partial_monomial_transpose(&final_mono_secret[i], &private_Q_inv[i]);
    }

    int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(mono_partial_secret[i].permutation[j]!=N+1){
                if(mono_partial_secret[i].coefficients[j]!=final_mono_secret[i].coefficients[j]
                    ||mono_partial_secret[i].permutation[j]!=final_mono_secret[i].permutation[j]){
                        fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                    res = failure;
             }
            }
        }
    }
    res = success;

    if(res == success)
        *num_cols += total_count_coeff;
    return res;

}

int SIM_recover_full_secret_key(int *num_sign){  
    
   pubkey_t pk;
   prikey_t sk;
   LESS_keygen(&sk,&pk);

   SHAKE_STATE_STRUCT sk_shake_state;
   initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
   unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
   for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
   }
   monomial_t final_mono_secret[NUM_KEYPAIRS-1];
   monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];
   
   for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
   }

   monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

   for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
   }

   int count_mono_coeff[NUM_KEYPAIRS-1]={0};
   int total_count_coeff = 0, count_iteration=0;

   while(total_count_coeff < N*(NUM_KEYPAIRS-1)){
      count_iteration++;
      sig_t signature;
      char message[8] = "Signme!";          /* We do not necessarily need same message */
      LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

      uint8_t fixed_weight_string[T] = {0};
      expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

      unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
      memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
      unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

      int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
      unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);
      
      monomial_t Q_tilde;
      normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
      #ifdef COMPRESS_CMT_COLUMNS
      uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
      #endif
      monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


      rref_generator_mat_t G0_rref;
      generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
      generator_mat_t full_G0;
      generator_rref_expand(&full_G0,&G0_rref);

      int employed_monoms = 0;
      monomial_action_IS_t mono_action;
      for(int i = 0; i < t; i++) {
         int F = fixed_weight_string[i];
         if ( F!= 0){
            if(count_mono_coeff[F-1]<N) {
               monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
               prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

               expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
               find_partial_monomial_update(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action);
               total_count_coeff-=count_mono_coeff[F-1];
               count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
               total_count_coeff+=count_mono_coeff[F-1];
            }  
            employed_monoms++;
         }
      }
   }

   //printf("The number of iteration is %d\n",count_iteration);
   for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_transpose(&final_mono_secret[i], &mono_partial_secret[i]);
   }
   
   int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(private_Q_inv[i].coefficients[j]!=final_mono_secret[i].coefficients[j]
             ||private_Q_inv[i].permutation[j]!=final_mono_secret[i].permutation[j]){
                    fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                res = failure;
             }
        }
    }
    res = success;

    if(res == success)
        *num_sign += count_iteration;
    return res;

}

int SIM_recover_full_secret_key_prob(int *num_sign){  
    
   pubkey_t pk;
   prikey_t sk;
   LESS_keygen(&sk,&pk);

   SHAKE_STATE_STRUCT sk_shake_state;
   initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
   unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
   for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
   }
   monomial_t final_mono_secret[NUM_KEYPAIRS-1];
   monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];
   
   for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
   }

   monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

   for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
   }

   int count_mono_coeff[NUM_KEYPAIRS-1]={0};
   int total_count_coeff = 0, count_iteration=0;
   int flag;

   while(total_count_coeff < N*(NUM_KEYPAIRS-1)){

      count_iteration++;

      flag = sample_binom_p();
      if(flag != 1)
        continue;

      sig_t signature;
      char message[8] = "Signme!";          /* We do not necessarily need same message */
      LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

      uint8_t fixed_weight_string[T] = {0};
      expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

      unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
      memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
      unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

      int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
      unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);
      
      monomial_t Q_tilde;
      normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
      #ifdef COMPRESS_CMT_COLUMNS
      uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
      #endif
      monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


      rref_generator_mat_t G0_rref;
      generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
      generator_mat_t full_G0;
      generator_rref_expand(&full_G0,&G0_rref);

      int employed_monoms = 0;
      monomial_action_IS_t mono_action;
      for(int i = 0; i < t; i++) {
         int F = fixed_weight_string[i];
         if ( F!= 0){
            if(count_mono_coeff[F-1]<N) {
               monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
               prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

               expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
               find_partial_monomial_update(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action);
               total_count_coeff-=count_mono_coeff[F-1];
               count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
               total_count_coeff+=count_mono_coeff[F-1];
            }  
            employed_monoms++;
         }
      }
   }

   //printf("The number of iteration is %d\n",count_iteration);
   for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_transpose(&final_mono_secret[i], &mono_partial_secret[i]);
   }
   
   int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(private_Q_inv[i].coefficients[j]!=final_mono_secret[i].coefficients[j]
             ||private_Q_inv[i].permutation[j]!=final_mono_secret[i].permutation[j]){
                    fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                res = failure;
             }
        }
    }
    res = success;

    if(res == success)
        *num_sign += count_iteration;
    return res;

}


/******************************************************
 * 
 * 
 * 
 *                   NEW SIMULATIONS 
 *              WITH IMPROVED KEY-RECOVERY
 * 
 * 
 * 
 * ****************************************************/

int SIM_secret_key_single_fault_new(int *num_cols){  
    
    pubkey_t pk;
    prikey_t sk;
    LESS_keygen(&sk,&pk);

    SHAKE_STATE_STRUCT sk_shake_state;
    initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
    unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
    for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
    }
    monomial_t final_mono_secret[NUM_KEYPAIRS-1];
    monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];

    for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
    }

    monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

    for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
    }

    int count_mono_coeff[NUM_KEYPAIRS-1]={0};
    int total_count_coeff = 0;

    sig_t signature;
    char message[8]; 
    csprng_randombytes(message,7,&sk_shake_state);
    message[8]='\0';

    LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

    uint8_t fixed_weight_string[T] = {0};
    expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

    unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
    memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
    unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

    int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
    unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

    monomial_t Q_tilde;
    normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
    #ifdef COMPRESS_CMT_COLUMNS
    uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
    #endif
    monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


    rref_generator_mat_t G0_rref;
    generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
    generator_mat_t full_G0;
    generator_rref_expand(&full_G0,&G0_rref);

    int employed_monoms = 0;
    monomial_action_IS_t mono_action;
    for(int i = 0; i < t; i++) {
        int F = fixed_weight_string[i];
        if ( F!= 0){
            //printf("%d\n",F);
            if(count_mono_coeff[F-1]<N) {
                monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
                prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

                expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
                find_partial_monomial_update(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action);

                // Get F-th public key

                generator_mat_t G_hat;
                expand_to_rref(&G_hat, pk.SF_G[F-1]);

                invertible_mat_t G_hat_J_star, G_0_J, S, G_tmp;

                int k=0;
                for(int x=0; x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            FQ_ELEM tmp = fq_inv(mono_partial_secret[F-1].coefficients[x]);
                            tmp = fq_red((FQ_DOUBLEPREC)tmp * (FQ_DOUBLEPREC)full_G0.values[j][x]);
                            G_0_J.values[j][k] = tmp;
                        }
                        k++;
                    }
                }

                k=0;
                for(int x=0;x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            G_hat_J_star.values[j][k] = G_hat.values[j][mono_partial_secret[F-1].permutation[x]];
                        }
                        k++;
                    }
                    
                }

                generator_mat_t G_left, G_right;

                G_tmp = G_0_J;

                matrix_inverse(&G_0_J);

                // Recovering S where G_hat = S * G_0 * Q^{-1}
                mul_invertible_invertible(&S, &G_hat_J_star, &G_0_J);
                matrix_inverse(&S);
                mul_invertible_generator(&G_left, &S, &G_hat);
                
                for(int x=0;x<N;x++)
                {
                    for(FQ_ELEM e = 1; e<Q; e++){
                        int j = compare(&G_left, &full_G0, e, x);
                        //if(j!=N+1) printf("%d %d\n", e, j);
                        if(j!=N+1) {
                            mono_partial_secret[F-1].permutation[x] = j;
                            mono_partial_secret[F-1].coefficients[x] = fq_inv(e);
                        }
                    }
                }

                count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
                total_count_coeff+=1;
                
            } 

            employed_monoms++;
        
        }
    }


    int res = success;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(mono_partial_secret[i].permutation[j]!=N+1){
                if(mono_partial_secret[i].coefficients[j]!=private_Q_inv[i].coefficients[j]
                    ||mono_partial_secret[i].permutation[j]!=private_Q_inv[i].permutation[j]){
                        //fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                    res = failure;
                }
            }
            else break;
        }
    }
    //res = success;

    if(res == success)
        *num_cols += total_count_coeff;
    return res;

}

int SIM_recover_full_secret_key_new(int *num_sign){  
    
   pubkey_t pk;
   prikey_t sk;
   LESS_keygen(&sk,&pk);

   SHAKE_STATE_STRUCT sk_shake_state;
   initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
   unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
   for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
   }
   monomial_t final_mono_secret[NUM_KEYPAIRS-1];
   monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];
   
   for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
   }

   monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

   for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
   }

   int count_mono_coeff[NUM_KEYPAIRS-1]={0};
   int total_count_coeff = 0, count_iteration=0;

   while(total_count_coeff < (NUM_KEYPAIRS-1)){
      count_iteration++;
      sig_t signature;
      char message[8] = "Signme!";          /* We do not necessarily need same message */
      LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

      uint8_t fixed_weight_string[T] = {0};
      expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

      unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
      memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
      unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

      int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
      unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);
      
      monomial_t Q_tilde;
      normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
      #ifdef COMPRESS_CMT_COLUMNS
      uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
      #endif
      monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


      rref_generator_mat_t G0_rref;
      generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
      generator_mat_t full_G0;
      generator_rref_expand(&full_G0,&G0_rref);

      int employed_monoms = 0;
      monomial_action_IS_t mono_action;
      for(int i = 0; i < t; i++) {
        int F = fixed_weight_string[i];
        if ( F!= 0){
            //printf("%d\n",F);
            if(count_mono_coeff[F-1]<N) {
                monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
                prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

                expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
                find_partial_monomial_update(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action);

                // Get F-th public key

                generator_mat_t G_hat;
                expand_to_rref(&G_hat, pk.SF_G[F-1]);

                invertible_mat_t G_hat_J_star, G_0_J, S, G_tmp;

                int k=0;
                for(int x=0; x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            FQ_ELEM tmp = fq_inv(mono_partial_secret[F-1].coefficients[x]);
                            tmp = fq_red((FQ_DOUBLEPREC)tmp * (FQ_DOUBLEPREC)full_G0.values[j][x]);
                            G_0_J.values[j][k] = tmp;
                        }
                        k++;
                    }
                }

                k=0;
                for(int x=0;x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            G_hat_J_star.values[j][k] = G_hat.values[j][mono_partial_secret[F-1].permutation[x]];
                        }
                        k++;
                    }
                    
                }

                generator_mat_t G_left, G_right;

                G_tmp = G_0_J;

                matrix_inverse(&G_0_J);

                // Recovering S where G_hat = S * G_0 * Q^{-1}
                mul_invertible_invertible(&S, &G_hat_J_star, &G_0_J);
                matrix_inverse(&S);
                mul_invertible_generator(&G_left, &S, &G_hat);
                
                for(int x=0;x<N;x++)
                {
                    for(FQ_ELEM e = 1; e<Q; e++){
                        int j = compare(&G_left, &full_G0, e, x);
                        //if(j!=N+1) printf("%d %d\n", e, j);
                        if(j!=N+1) {
                            mono_partial_secret[F-1].permutation[x] = j;
                            mono_partial_secret[F-1].coefficients[x] = fq_inv(e);
                        }
                    }
                }

                
                count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
                total_count_coeff+=1;
                
            } 

            employed_monoms++;
        
        }
      }
   }

   //printf("The number of iteration is %d\n",count_iteration);
   
   int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(private_Q_inv[i].coefficients[j]!=mono_partial_secret[i].coefficients[j]
             ||private_Q_inv[i].permutation[j]!=mono_partial_secret[i].permutation[j]){
                    fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                res = failure;
             }
        }
    }
    res = success;

    if(res == success)
        *num_sign += count_iteration;
    return res;

}

int SIM_recover_full_secret_key_prob_new(int *num_sign){  
    
   pubkey_t pk;
   prikey_t sk;
   LESS_keygen(&sk,&pk);

   SHAKE_STATE_STRUCT sk_shake_state;
   initialize_csprng(&sk_shake_state,sk.compressed_sk,SEED_LENGTH_BYTES);

     /* generation, hence NUM_KEYPAIRS-1 */
   unsigned char monomial_seeds[NUM_KEYPAIRS-1][SEED_LENGTH_BYTES];
   for (int i = 0; i < NUM_KEYPAIRS-1; i++) {
      csprng_randombytes(monomial_seeds[i],
                           SEED_LENGTH_BYTES,
                           &sk_shake_state);
   }
   monomial_t final_mono_secret[NUM_KEYPAIRS-1];
   monomial_t private_Q_transpose[NUM_KEYPAIRS-1], private_Q_inv[NUM_KEYPAIRS-1];
   
   for(int i = 0; i < NUM_KEYPAIRS-1; i++) {
      monomial_mat_seed_expand(&private_Q_inv[i], monomial_seeds[i]);
   }

   monomial_t mono_partial_secret[NUM_KEYPAIRS-1];

   for(int i=0; i<NUM_KEYPAIRS-1; i++){
      partial_monomial_initialize(&final_mono_secret[i]);
      partial_monomial_initialize(&mono_partial_secret[i]);
   }

   int count_mono_coeff[NUM_KEYPAIRS-1]={0};
   int total_count_coeff = 0, count_iteration=0;
   int flag;

   while(total_count_coeff < (NUM_KEYPAIRS-1)){

      count_iteration++;

      flag = sample_binom_p();
      if(flag != 1)
        continue;

      sig_t signature;
      char message[8] = "Signme!";          /* We do not necessarily need same message */
      LESS_sign(&sk,message,8,&signature);  /* Fault is induced at line 102-104 of Reference_implementation/lib/seedtree.c */

      uint8_t fixed_weight_string[T] = {0};
      expand_digest_to_fixed_weight(fixed_weight_string,signature.digest);

      unsigned char ephem_monomials_seed[SEED_LENGTH_BYTES];
      memcpy(ephem_monomials_seed, signature.seed_storage, SEED_LENGTH_BYTES);
      unsigned char seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES] = {0};

      int t = generate_partial_seed_tree( seed_tree, ephem_monomials_seed,signature.tree_salt);
      unsigned char *ephem_monomial_seeds = seed_tree + SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);
      
      monomial_t Q_tilde;
      normalized_IS_t* V_array = calloc(T, sizeof(normalized_IS_t));
      #ifdef COMPRESS_CMT_COLUMNS
      uint8_t *V_array_compressed = calloc(T*RREF_IS_COLUMNS_PACKEDBYTES, sizeof(uint8_t));
      #endif
      monomial_action_IS_t* Q_bar_actions = calloc(T, sizeof(monomial_action_IS_t));


      rref_generator_mat_t G0_rref;
      generator_SF_seed_expand(&G0_rref, pk.G_0_seed);
      generator_mat_t full_G0;
      generator_rref_expand(&full_G0,&G0_rref);

      int employed_monoms = 0;
      monomial_action_IS_t mono_action;
      for(int i = 0; i < t; i++) {
         int F = fixed_weight_string[i];
         if ( F!= 0){
            //printf("%d\n",F);
            if(count_mono_coeff[F-1]<N) {
                monomial_mat_seed_expand(&Q_tilde, ephem_monomial_seeds+i*SEED_LENGTH_BYTES);
                prepare_digest_input(&V_array[i],&Q_bar_actions[i],&full_G0, &Q_tilde);

                expand_to_monom_action(&mono_action, signature.monom_actions[employed_monoms]);
                find_partial_monomial_update(&mono_partial_secret[F-1], &Q_bar_actions[i], &mono_action);

                // Get F-th public key

                generator_mat_t G_hat;
                expand_to_rref(&G_hat, pk.SF_G[F-1]);

                invertible_mat_t G_hat_J_star, G_0_J, S, G_tmp;

                int k=0;
                for(int x=0; x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            FQ_ELEM tmp = fq_inv(mono_partial_secret[F-1].coefficients[x]);
                            tmp = fq_red((FQ_DOUBLEPREC)tmp * (FQ_DOUBLEPREC)full_G0.values[j][x]);
                            G_0_J.values[j][k] = tmp;
                        }
                        k++;
                    }
                }

                k=0;
                for(int x=0;x<N;x++){
                    if(mono_partial_secret[F-1].permutation[x] != N+1){
                        for(int j=0;j<K;j++){
                            G_hat_J_star.values[j][k] = G_hat.values[j][mono_partial_secret[F-1].permutation[x]];
                        }
                        k++;
                    }
                    
                }

                generator_mat_t G_left, G_right;

                G_tmp = G_0_J;

                matrix_inverse(&G_0_J);

                // Recovering S where G_hat = S * G_0 * Q^{-1}
                mul_invertible_invertible(&S, &G_hat_J_star, &G_0_J);
                matrix_inverse(&S);
                mul_invertible_generator(&G_left, &S, &G_hat);
                
                for(int x=0;x<N;x++)
                {
                    for(FQ_ELEM e = 1; e<Q; e++){
                        int j = compare(&G_left, &full_G0, e, x);
                        //if(j!=N+1) printf("%d %d\n", e, j);
                        if(j!=N+1) {
                            mono_partial_secret[F-1].permutation[x] = j;
                            mono_partial_secret[F-1].coefficients[x] = fq_inv(e);
                        }
                    }
                }

                count_mono_coeff[F-1]=update_count_mono(&mono_partial_secret[F-1]);
                total_count_coeff+=1;
                
            } 

            employed_monoms++;
        
        }
      }
   }

   //printf("The number of iteration is %d\n",count_iteration);
   
   int res;
   /* Checking if the recovered secret is same as the actual secret */
    for(int i=0; i<NUM_KEYPAIRS-1; i++){
        for(int j=0; j<N; j++){
            if(private_Q_inv[i].coefficients[j]!=mono_partial_secret[i].coefficients[j]
             ||private_Q_inv[i].permutation[j]!=mono_partial_secret[i].permutation[j]){
                    fprintf(stderr,"Recovered secret-key is not same as original secret-key\n");
                res = failure;
             }
        }
    }
    res = success;

    if(res == success)
        *num_sign += count_iteration;
    return res;

}
