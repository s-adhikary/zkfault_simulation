/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "timing_and_stat.h"
#include "fq_arith.h"
#include "seedtree.h"
//#include "merkle_tree.h"
#include <assert.h>
#include "pack_unpack.h"
#define NUM_TEST_ITERATIONS 100
#include "restr_arith.h"
#include "arith_unit_tests.h"
#include "CROSS.h"
#include "csprng_hash.h"
#include "parameters.h"

#define LEFT_CHILD(i) (2*i+1)
#define RIGHT_CHILD(i) (2*i+2)
#define PARENT(i) ((i-1)/2)
#define SIBLING(i) ( ((i)%2) ? i+1 : i-1 )
#define IS_LEFT_SIBLING(i) (i%2)    
//#define RSDP 1

void info(void){
    fprintf(stderr,"CROSS functional testing utility\n");
#if defined(RSDP)
    fprintf(stderr,"RSDP Variant\n");
#elif defined(RSDPG)
    fprintf(stderr,"RSDPG Variant\n");
#endif
    fprintf(stderr,"Code parameters: n= %d, k= %d, q=%d\n", N,K,Q);
    fprintf(stderr,"restriction size: z=%d\n",Z);
    fprintf(stderr,"Fixed weight challenge vector: %d rounds, weight %d \n",T,W);
    fprintf(stderr,"Private key: %luB\n", sizeof(prikey_t));
    fprintf(stderr,"Public key %luB\n", sizeof(pubkey_t));
    fprintf(stderr,"Signature: %luB\n", sizeof(sig_t));
}

//#endif


static
void expand_public_seed1(FQ_ELEM V_tr[N-K][K],
                        const uint8_t seed_pub[SEED_LENGTH_BYTES]){
  CSPRNG_STATE_T CSPRNG_state_mat;
  initialize_csprng(&CSPRNG_state_mat, seed_pub, SEED_LENGTH_BYTES);
  CSPRNG_fq_mat(V_tr,&CSPRNG_state_mat);
}

static
void expand_private_seed1(FZ_ELEM eta[N], FQ_ELEM V_tr[N-K][K], const uint8_t seed[SEED_LENGTH_BYTES]){
  uint8_t seede_seed_pub[2][SEED_LENGTH_BYTES];
  CSPRNG_STATE_T CSPRNG_state;
  initialize_csprng(&CSPRNG_state,seed,SEED_LENGTH_BYTES);
  csprng_randombytes((uint8_t *)seede_seed_pub,
                     2*SEED_LENGTH_BYTES,
                     &CSPRNG_state);

  expand_public_seed1(V_tr,seede_seed_pub[1]);

  CSPRNG_STATE_T CSPRNG_state_eta;
  initialize_csprng(&CSPRNG_state_eta, seede_seed_pub[0], SEED_LENGTH_BYTES);
  CSPRNG_zz_vec(eta,&CSPRNG_state_eta);
}


/* returns 1 if the test is successful, 0 otherwise */
int CROSS_sign_verify_test(){
    pubkey_t pk;
    prikey_t sk;
    sig_t signature;
    char message[8] = "Signme!";
    CROSS_keygen(&sk,&pk);
    CROSS_sign(&sk,message,8,&signature);
    int is_signature_ok;
    is_signature_ok = CROSS_verify(&pk,message,8,&signature);
    return is_signature_ok;
}
//new function
int regenerate_leaves_update(unsigned char
                      seed_tree[NUM_NODES_OF_SEED_TREE*SEED_LENGTH_BYTES],
                      const unsigned char indices_to_publish[T],
                      const unsigned char *stored_seeds,
                      const unsigned char salt[SALT_LENGTH_BYTES])
{
   /* complete linearized binary tree containing boolean values determining
    * if a node is to be released or not. Nodes set to 1 are not to be released
    * oldest ancestor of sets of nodes equal to 0 are to be released */
   //unsigned char flags_tree_to_publish[2*NUM_LEAVES_OF_SEED_TREE-1] = {0};
   //compute_seeds_to_publish(flags_tree_to_publish, indices_to_publish);

   const uint32_t csprng_input_len = SALT_LENGTH_BYTES +
                                     SEED_LENGTH_BYTES +
                                     sizeof(uint32_t);
   unsigned char csprng_input[csprng_input_len];
   CSPRNG_STATE_T tree_csprng_state;

  //  memcpy(csprng_input, salt, HASH_DIGEST_LENGTH);

   memcpy(csprng_input, salt, SALT_LENGTH_BYTES);
   uint32_t k=1,count=0;
   memcpy(seed_tree+SEED_LENGTH_BYTES, stored_seeds,SEED_LENGTH_BYTES);
   /*printf("Seed 1: \n");
   for(int j=0;j<SEED_LENGTH_BYTES; j++)
    printf("%hx ", (seed_tree+SEED_LENGTH_BYTES)[j]);
    printf("\n");*/
    memcpy(csprng_input + SEED_LENGTH_BYTES,
                   seed_tree + SEED_LENGTH_BYTES,
                   SEED_LENGTH_BYTES);
       //  printf("%d\n", i);   
      *((uint32_t *)(csprng_input + SEED_LENGTH_BYTES + SEED_LENGTH_BYTES)) = 1;
      initialize_csprng(&tree_csprng_state, csprng_input, csprng_input_len);
      csprng_randombytes(seed_tree + LEFT_CHILD(1)*SEED_LENGTH_BYTES,
                         2*SEED_LENGTH_BYTES,
                         &tree_csprng_state);

   //int nodes_used = 0;
   /* regenerating the seed tree never starts from the root, as it is never
    * disclosed*/
   //printf("NUM_LEAVES_OF_SEED_TREE: %d\n", NUM_LEAVES_OF_SEED_TREE);
   for (uint32_t i = 3; i < NUM_LEAVES_OF_SEED_TREE-1; ) {
      /* if the current node is a seed which was published, memcpy it in place */
      memcpy(csprng_input + SEED_LENGTH_BYTES,
                   seed_tree + i*SEED_LENGTH_BYTES,
                   SEED_LENGTH_BYTES);
       // printf("%d\n", i);   
      *((uint32_t *)(csprng_input + SEED_LENGTH_BYTES + SEED_LENGTH_BYTES)) = i;
      initialize_csprng(&tree_csprng_state, csprng_input, csprng_input_len);
      csprng_randombytes(seed_tree + LEFT_CHILD(i)*SEED_LENGTH_BYTES,
                         2*SEED_LENGTH_BYTES,
                         &tree_csprng_state);
        /*printf("Seed %d: \n", 2*i);
   for(int j=0;j<2*SEED_LENGTH_BYTES; j++)
        printf("%hx ", (seed_tree + (LEFT_CHILD(i)*SEED_LENGTH_BYTES))[j]);
    printf("\n");*/

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
} /* end regenerate_leaves */



void Attack_CROSS(){
    int t=0;
    pubkey_t pk;
    prikey_t sk;
    sig_t signature;
    char message[8] = "Signme!";
     
    // Key generation of victim device
    CROSS_keygen(&sk,&pk);
    // that component is common in pk and sk, so we can ignore this
    FQ_ELEM V_tr1[N-K][K]; 
    FQ_ELEM V_tr2[N-K][K];
    // This is the secret component that we want to find and later we will use this for prove of correctness
    FZ_ELEM eta1[N]; 
    FZ_ELEM eta2[N]; 
    #if defined(RSDP)
        expand_private_seed1(eta2,V_tr2,sk.seed);
    #elif defined(RSDPG)
        FZ_ELEM zeta[M];
    FZ_ELEM W_tr[N-M][M];
    expand_private_seed(eta1,zeta,V_tr1,W_tr,sk.seed);
    #endif

    // Call the Faulted CROSS_sign Algorithm
    CROSS_sign(&sk,message,8,&signature);


    // The method of finding secret component eta1
    uint8_t fixed_weight_b[T]={0};
    expand_digest_to_fixed_weight(fixed_weight_b,signature.digest_b);
     for(int i = 0; i<T; i++){
        fixed_weight_b[i] = !fixed_weight_b[i];
     }

    uint8_t seed_tree[SEED_LENGTH_BYTES*NUM_NODES_OF_SEED_TREE]={0};
    regenerate_leaves_update(seed_tree,fixed_weight_b, signature.stp,signature.salt);
    uint8_t * seed_tree_leaves = seed_tree +
                                 SEED_LENGTH_BYTES*(NUM_LEAVES_OF_SEED_TREE-1);

    uint8_t cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES+sizeof(uint16_t)];
    memcpy(cmt_1_i_input+SEED_LENGTH_BYTES, signature.salt, SALT_LENGTH_BYTES);
    uint8_t cmt_1[T][HASH_DIGEST_LENGTH] = {0};

    uint8_t beta_buf[2*HASH_DIGEST_LENGTH+SALT_LENGTH_BYTES];
    hash(beta_buf,(uint8_t*) message,8);
    memcpy(beta_buf+HASH_DIGEST_LENGTH, signature.digest_01, HASH_DIGEST_LENGTH);
    memcpy(beta_buf+2*HASH_DIGEST_LENGTH, signature.salt, SALT_LENGTH_BYTES);
     uint8_t d_beta[HASH_DIGEST_LENGTH];
    hash(d_beta,beta_buf,sizeof(beta_buf));
    CSPRNG_STATE_T CSPRNG_state;
    FQ_ELEM beta[T];
    initialize_csprng(&CSPRNG_state,d_beta,HASH_DIGEST_LENGTH);
    CSPRNG_fq_vec_beta(beta, &CSPRNG_state);
    FZ_ELEM eta_tilde[N];
    #if defined(RSDP)
    uint8_t cmt_0_i_input[sizeof(FQ_ELEM)*(N-K)+
                          sizeof(FZ_ELEM)*N+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];
    const int offset_salt = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*N;
//    const int offset_round_idx = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*N+SALT_LENGTH_BYTES;
#elif defined(RSDPG)
printf("Hii\n");
  //  FZ_ELEM zeta_tilde[M];
   // FZ_ELEM delta[T][M];
    uint8_t cmt_0_i_input[sizeof(FQ_ELEM)*(N-K)+
                          sizeof(FZ_ELEM)*M+
                          SALT_LENGTH_BYTES+sizeof(uint16_t)];
    const int offset_salt = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*M;
  //  const int offset_round_idx = sizeof(FQ_ELEM)*(N-K)+sizeof(FZ_ELEM)*M+SALT_LENGTH_BYTES;
#endif

    memcpy(cmt_0_i_input+offset_salt, signature.salt, SALT_LENGTH_BYTES);
    FQ_ELEM y[T][N];
    int used_rsps = 0;
    uint8_t seed_u_t_seed_e_t[2*SEED_LENGTH_BYTES];
  //  int is_signature_ok = 1;
    //printf("Attack funcion\n");
    for(int i = 0; i<(NUM_LEAVES_OF_SEED_TREE>>1); i++){
        // We run the if loop only once.
        if(fixed_weight_b[i] != 1)
        {
          //  Compute eta_tilde
            memcpy(cmt_1_i_input,
                   seed_tree_leaves+SEED_LENGTH_BYTES*i,
                   SEED_LENGTH_BYTES);
            cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES] = (i >> 8) &0xFF;
            cmt_1_i_input[SEED_LENGTH_BYTES+SALT_LENGTH_BYTES+1] = i & 0xFF; 
            hash(cmt_1[i],cmt_1_i_input,sizeof(cmt_1_i_input));
            initialize_csprng(&CSPRNG_state,
                              seed_tree_leaves+SEED_LENGTH_BYTES*i,
                              SEED_LENGTH_BYTES);
            csprng_randombytes(seed_u_t_seed_e_t,
                               2*SEED_LENGTH_BYTES,
                               &CSPRNG_state);
            #if defined(RSDP)
            /* expand eta_tilde */
            initialize_csprng(&CSPRNG_state,
                              seed_u_t_seed_e_t+SEED_LENGTH_BYTES,
                              SEED_LENGTH_BYTES);
            CSPRNG_zz_vec(eta_tilde, &CSPRNG_state);
            
            #elif defined(RSDPG)
            initialize_csprng(&CSPRNG_state,
                              seed_u_t_seed_e_t+SEED_LENGTH_BYTES,
                              SEED_LENGTH_BYTES);
            FZ_ELEM zeta_tilde[M];
            CSPRNG_zz_inf_w(zeta_tilde, &CSPRNG_state);
            fz_inf_w_by_fz_matrix(eta_tilde,zeta_tilde,W_tr);
            fz_dz_norm_sigma(eta_tilde);
            printf("Eta_tilde in attack %d\n",i);
            for(int j=0; j<N; j++)\
                printf("%hd ", eta_tilde[j]);
            printf("\n");
            #endif 
            

// Compute sigma
            unpack_fq_vec(y[i], signature.rsp_0[used_rsps].y);
            FZ_ELEM sigma_local[N];
            
            #if defined(RSDP)
            /*sigma is memcpy'ed directly into cmt_0 input buffer */
            FZ_ELEM* sigma_ptr = cmt_0_i_input+sizeof(FQ_ELEM)*(N-K);
	        unpack_fz_vec(sigma_local, signature.rsp_0[used_rsps].sigma);
            #elif defined(RSDPG)
            int is_signature_ok = 1;
            /*delta is memcpy'ed directly into cmt_0 input buffer */
            FZ_ELEM* sigma_ptr = cmt_0_i_input+sizeof(FQ_ELEM)*(N-K);
	        unpack_fz_rsdp_g_vec(sigma_ptr, signature.rsp_0[used_rsps].delta);
            is_signature_ok = is_signature_ok &&
                              is_zz_inf_w_valid(sigma_ptr);
            fz_inf_w_by_fz_matrix(sigma_local,sigma_ptr,W_tr);
            printf("Sigma_local in attack function %d\n",i);
            for(int j=0; j<N; j++)\
            printf("%hd ", sigma_local[j]);
        printf("\n");
        #endif
            
            used_rsps++;

            //Compute the secret component eta = sigma_local+eta_tilde mod 2^{3}-1 
            FZ_ELEM eta[N];
            int count=0;
            #if defined(RSDP) 
            for(int j = 0; j < N; j++){
              eta[j]= (( sigma_local[j] + eta_tilde[j] )&0x07)+(( sigma_local[j] + eta_tilde[j] )>>3);
               if(eta[j]==7)
                eta[j]=0;
            }
            for(int j=0; j<N; j++)
            {
                if(eta[j]!=eta2[j]){
                    printf("Print the value does not match for %d th location: eta[%d]:%hd, eta1[%d]:%hd \n", j, j, eta[j], eta2[j]);
                    count++;
                }
            }
            #elif defined(RSDPG)
            for(int j = 0; j < N; j++){
              eta[j]= (( sigma_local[j] + eta_tilde[j] )&0x7f)+(( sigma_local[j] + eta_tilde[j] )>>7);
               if(eta[j]==127)
                eta[j]=0;
            }
             for(int j=0; j<N; j++)
            {
                if(eta[j]!=eta1[j]){
                    printf("Print the value does not match for %d th location: eta[%d]:%hd, eta1[%d]:%hd \n", j, j, eta[j], eta1[j]);
                    count++;
                }
            }
            #endif
            // Check whether eta =? eta1 or not
            
            #if defined(RSDP)
            if(count==0)
                printf("The attack is successful for RSDP\n");
            #elif defined(RSDPG)
            if(count==0)
                printf("The attack is successful for RSDPG\n");
            #endif
            break;
        }
    }
  //  return t;
}
int main(int argc, char* argv[]){
    initialize_csprng(&platform_csprng_state,
                      (const unsigned char *)"012345678912345",
                      16);
    fprintf(stderr,"CROSS reference implementation functional testbench\n");
    info();
    int tests_ok = 0;
    for (int i = 0; i < NUM_TEST_ITERATIONS; i++) {
            Attack_CROSS();       
    }
    fprintf(stderr,"\n%d tests functional out of %d\n",tests_ok,NUM_TEST_ITERATIONS);
    return 0;
}
