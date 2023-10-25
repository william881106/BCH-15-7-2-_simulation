#include <iostream>
#include <vector>
#include <random>
#include <cstdlib> /* 亂數相關函數 */
#include <ctime>   /* 時間相關函數 */

#define coding_len 15
#define info_len 7
#define error_correct 2
#define minimun_poly_order 8
using namespace std;

#pragma region initial_setting
//-------------------------------------------------------------------------------------------------
//g(x) = 1 + x^4 + x^6 + x^7 + x^8  ---------------------------------------------------------------
//p(x) = 1 + x   + x^4
//-------------------------------------------------------------------------------------------------

int prime_poly = 19;  //p(x) = 1 + x + x^4 , 10011的十進位
int alpha_reg = 1;    

vector<int> generator_poly = {1,0,0,0,1,0,1,1,1};
vector< vector<int> > GF_table_binary(16, vector<int>(4,0));
vector<int> GF_table_decimal(16,0);
vector< vector<int> > G(7, vector<int>(15,0));


/* 標準常態分布 */
random_device rd;
mt19937 generator( rd() );
normal_distribution<double> norm(0.0,1.0);


int trials = 100000;
int error = 0;
int error_min = 100;
double coderate = (1.0*info_len)/(1.0*coding_len);
double randn_num;
double No;
double noise_temp;

vector<int> message(7,0);
vector<int> codeword(15,0);
vector<double> noise(15,0);
vector<double> y(15,0);
vector<int> y_hat(15,0);


int find_symdrome(vector<int> y_hat, int alpha_power, vector<int> G);
vector<int> symdrome_decimal(2*error_correct,0);
vector<int> symdrome_power(2*error_correct,0);


vector<int> fifteen_vector_for_check(symdrome_decimal.size(), 15);
vector<double> u(error_correct+2,0);
vector<vector<int>> sigma_x(error_correct+2,vector<int>(minimun_poly_order,15));
vector<int> du(error_correct+2,0);
vector<int> lu(error_correct+2,0);
vector<int> u2_lu(error_correct+2,0);
int find_next_du(int current_step, vector<int> sigma_x, vector<int> symdrome_power, vector<int> GF_table_decimal, int lu);


int find_rho(int current_step, vector<int> u2_lu, vector<double> u, vector<int> du);
vector<int> find_next_sigma_x(int current_step, int rho, vector<vector<int>> sigma_x, vector<double> u, vector<int> du, vector<int> symdrome_power, vector<int> GF_table_decimal);
vector<int> final_sigma_x_power(5,0);
vector<int> error_patten(5,0);


vector<int> decode_polynomial(7,0);

#pragma endregion initial_setting

int main()
{   
    srand( time (0) );
    #pragma region construct_Galois_Field
    //-------------------------------------------------------------------------------------------------
    //construction of Galois Field---------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------
    GF_table_binary[15][0] = 0;
    GF_table_binary[15][1] = 0;
    GF_table_binary[15][2] = 0;
    GF_table_binary[15][3] = 0;

    for(int i = 0; i < coding_len; i++)
    {
        //此表格為0、1、a、a^2... mod(prime_poly) 的結果
        //prime polynomial order = 4，則取mod的order不會 >= 4
        //也就是遇到 10000 (16)要做判斷

        if( (alpha_reg & coding_len+1) == coding_len+1 )
        {
            alpha_reg = alpha_reg ^ prime_poly;
            //printf("%d\n", alpha_reg ^ prime_poly);
        }
        GF_table_binary[i][0] =  alpha_reg     & 1;
        GF_table_binary[i][1] = (alpha_reg>>1) & 1;
        GF_table_binary[i][2] = (alpha_reg>>2) & 1;
        GF_table_binary[i][3] = (alpha_reg>>3) & 1;
        
        alpha_reg = alpha_reg << 1;
        //printf("%d\n", alpha_reg);
    }

    for (int i = 0; i <= coding_len; i++) //改成十進位方便後續計算
    {
        for (int j = 0; j < 3; j++)
        {
            GF_table_decimal[i] = GF_table_binary[i][0] + 2*GF_table_binary[i][1] + 4*GF_table_binary[i][2] + 8*GF_table_binary[i][3];
        }
    }

    /*
    for(int i = 0; i <= coding_len; i++)
    {
        printf("%d%d%d%d\n", GF_table_binary[i][0], GF_table_binary[i][1], GF_table_binary[i][2], GF_table_binary[i][3]);
        cout << GF_table_decimal[i] << endl ;
    }
    cout<< endl;
    */

    #pragma endregion construct_Galois_Field

    #pragma region construct_G
    //-------------------------------------------------------------------------------------------------
    //construct G matrix (size k*n = 7*15)-------------------------------------------------------------
    //using Guass elimination to make G = [P|I] (systematic form)--------------------------------------
    //using systematic form for easier decoding--------------------------------------------------------
    //-------------------------------------------------------------------------------------------------

    for(int i = 0; i < info_len; i++) //construct G matrix
    {
        for(int j = 0; j <= coding_len-info_len; j++)
            G[i][j+i] = generator_poly[j];

    }
    
    for (int i = 0; i < info_len; i++) //Guass elimination G to Gsys
    {
        for (int ii = i+1; ii < info_len; ii++)
        {
            if(G[ii][coding_len-info_len + i])
            {
                for (int k = 0; k < coding_len; k++)
                {
                    G[ii][k] = G[i][k] ^ G[ii][k];
                }
            }
            else
                continue;
        }
    }
    
    /*   
    for(int i = 0; i < info_len; i++)
    {
        for(int j = 0; j < coding_len; j++)
            cout << G[i][j] <<" ";
        cout << "\n";
    }
    */ 

    #pragma endregion construct_G

    //-------------------------------------------------------------------------------------------------
    //simulation： SNR 0~12, --------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------
    
    for (int SNR_dB = 0; SNR_dB <= 24; SNR_dB++)
    {   
        double SNR_dB_2 = SNR_dB / 2.0;
        No = pow( 10.0, double(-SNR_dB_2) / 10.0 );
        error = 0;
        int times = 0;
        //for (int qq = 0; qq < trials; qq++)  
        while(times < trials || error < error_min)//||error_count<error_min	
        {
            times += 1;

            #pragma region generate_random_7_bits_data
            //-------------------------------------------------------------------------------------------------
            //generate random 7-bits data ---------------------------------------------------------------------
            //-------------------------------------------------------------------------------------------------
            fill(message.begin(), message.end(), 0);
            for(int i=0; i < info_len; i++)  
                message[i] = rand() %2;

            /*
            for(int i=0; i < info_len; i++) 
                cout << message[i] << " ";
            cout << endl;
            */
            #pragma endregion generate_random_7_bits_data

            #pragma region encode_data
            //-------------------------------------------------------------------------------------------------
            //encode  data ------------------------------------------------------------------------------------
            //[m0 m1 m2 m3 m4 m5 m6 m7] * G -------------------------------------------------------------------
            //-------------------------------------------------------------------------------------------------
            fill(codeword.begin(), codeword.end(), 0);
            for (int j = 0; j < coding_len; j++)
            {
                for (int i = 0; i < info_len; i++)
                {
                    codeword[j] += message[i]*G[i][j];         
                }
                codeword[j] = codeword[j] % 2;
            }
            /*
            cout << "encode data:       ";
            for(int i=0; i < coding_len; i++)
                cout << codeword[i] << " ";
            cout << endl;
            */
            #pragma endregion encode_data

            #pragma region generate_noise
            //-------------------------------------------------------------------------------------------------
            //generate noise  ---------------------------------------------------------------------------------
            //-------------------------------------------------------------------------------------------------
            for (int j = 0; j < coding_len; j++)
            {
                randn_num = norm(generator);
                noise[j] = randn_num * sqrt(No/2) * sqrt(1.0/coderate) ;
                //cout << SNR_dB << " : " << No<< " - " << noise[j] << endl;
            }
            #pragma endregion generate_noise

            #pragma region BPSK_modulation
            //-------------------------------------------------------------------------------------------------
            //BPSK_modulation with 1 / -1 ---------------------------------------------------------------------
            //receive signal y = m + n ------------------------------------------------------------------------
            //-------------------------------------------------------------------------------------------------
            for (int j = 0; j < coding_len; j++)
            {
                y[j] = ( 2.0 * codeword[j] - 1 ) + noise[j];
                //cout << SNR_dB << " : " << No << " " << y[j] << endl;
            }
            #pragma endregion BPSK_modulation

            #pragma region hard_decision

            //-------------------------------------------------------------------------------------------------
            //hard decision -----------------------------------------------------------------------------------
            //-------------------------------------------------------------------------------------------------
            //cout << "hard decision data:";
            for (int j = 0; j < coding_len; j++)
            {
                if( y[j] >= 0 ) 
                    y_hat[j] = 1; 
                else
                    y_hat[j] = 0; 
                //cout <<  y_hat[j] << " ";
            }
            //cout << endl;

            /*
            y_hat = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1 };
            for (int j = 0; j < coding_len; j++)
            {
                cout <<  y_hat[j] << " ";
            }
            cout << endl;
            */

            #pragma endregion hard_decision

            #pragma region find_symdrome

            //-------------------------------------------------------------------------------------------------
            //find symdrome -----------------------------------------------------------------------------------
            //-------------------------------------------------------------------------------------------------
            fill(symdrome_decimal.begin(), symdrome_decimal.end(), 0);
            fill(symdrome_power.begin(), symdrome_power.end(), 0);

            for (int i = 0; i < 2*error_correct; i++)
            {
                symdrome_decimal[i] = find_symdrome(y_hat, i+1, GF_table_decimal);
                //cout << symdrome_decimal[i] <<" " ;
            }
            //cout << endl;
            for (int i = 0; i < 2*error_correct; i++)
            {   
                for (int j = 0; j <= coding_len; j++)
                {
                    if (symdrome_decimal[i] == GF_table_decimal[j])
                        symdrome_power[i] = j;
                }
                //cout << symdrome_power[i] <<" " ;
            }
            //cout << endl;

            #pragma endregion find_symdrome

            
            //-------------------------------------------------------------------------------------------------
            //Berlekamps iterative algorithm ------------------------------------------------------------------
            //if all symdrome = 0 , means no error, than skip this algorithm ----------------------------------
            //-------------------------------------------------------------------------------------------------
            
            if( symdrome_power ==  fifteen_vector_for_check)
            {
                //cout << "no error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl;
            }
            else
            {   
                #pragma region Berlekamps
                //---------------initial parameter of the table -----------------------------------------------
                u[0] = -0.5;
                u[1] = 0;
                u[2] = 1;
                u[3] = 2;

                sigma_x[0][0] = 0; //power of alpha, not means number"0" :   if power = 15 means real number"0"
                sigma_x[1][0] = 0; //power of alpha

                du[0] = 0; //power of alpha
                du[1] = symdrome_power[0]; 

                lu[0] = 0;
                lu[1] = 0; 

                u2_lu[0] = -1;
                u2_lu[1] = 0; 

                int rho = 0;
                //---------------Berlekamps-----------------------------------------------
                
                for (int step = 1; step < error_correct; step++)
                {
                    if ( du[step] == 15 )
                    {
                        sigma_x[step+1] = sigma_x[step];
                        lu[step+1] = lu[step];       
                        du[step+1] = find_next_du(step, sigma_x[step+1], symdrome_power, GF_table_decimal,lu[step+1]);
                        u2_lu[step+1] = 2*u[step+1] - lu[step+1];
                    }
                    else
                    {
                        rho = find_rho(step, u2_lu, u, du);
                        sigma_x[step+1] = find_next_sigma_x(step, rho, sigma_x, u, du, symdrome_power, GF_table_decimal);
                        for (int i = sigma_x[step+1].size()-1; i > 0; i--)
                        {
                            if (sigma_x[step+1][i]==15)
                            {

                            }
                            else
                                lu[step+1] = i;
                        }
                        
                        du[step+1] = find_next_du(step, sigma_x[step+1], symdrome_power, GF_table_decimal, lu[step+1]);
                        u2_lu[step+1] = 2*u[step+1] - lu[step+1];
                        
                    }
                }
                
                //---------------Berlekamps Final update-----------------------------------------------
                int step = error_correct;
                if ( du[step] == 15 )
                    sigma_x[step+1] = sigma_x[step];
                else
                {
                    rho = find_rho(step, u2_lu, u, du);
                    sigma_x[step+1] = find_next_sigma_x(step, rho, sigma_x, u, du, symdrome_power, GF_table_decimal);
                }

                final_sigma_x_power = sigma_x[step+1];


                


                #pragma endregion Berlekamps

                #pragma region find_error_position
            
                int t = 0;
                int power = 0;
                int s = 0;
                for(int i = 0; i < coding_len; i++)//從a^0試到a^14
                {

                    s = GF_table_decimal[final_sigma_x_power[0]];
                    //cout << s;
                    for(int j = 1; j <= error_correct; j++)
                    {
                        if (final_sigma_x_power[j] != 15)
                        {
                            power = (final_sigma_x_power[j] + i*j) % coding_len;
                            s = s ^ GF_table_decimal[power];
                        }
                        else
                            continue;
                    }

                    if(s == 0)
                    {
                        //cout << endl << i <<endl;
                        error_patten[t] = (coding_len-i) % coding_len;//計算錯誤的點在哪裡 
                        
                        t = t + 1;
                    }
                    
                    if(t==error_correct)
                    {
                        
                        break;
                    }
                }
                //cout << endl <<  error_patten[0] <<" "<<  error_patten[1] <<endl;

                #pragma endregion find_error_position

                #pragma region correct_error

                for(int i = 0; i < t; i++)
                {
                    y_hat[error_patten[i]] = y_hat[error_patten[i]]^1;
                }    
                
                //cout << "decode data:       ";
                for(int i = 0; i < coding_len; i++)  
                {  
                    //cout <<  y_hat[i] << " ";
                }
                //cout << endl;    
                

                #pragma endregion correct_error
                
                #pragma region cout_for_checking_Berlekamps
                
                cout << "-----------------------------------------------------------------------------"<< endl;
                cout << " ---u--- "<<" --------sigma_x-------- "<<" ---du--- "<<" ---lu--- "<<" ---2u-lu--- "<<" ---------"<< endl;
                cout << "-----------------------------------------------------------------------------"<< endl;
                cout << "  "<<u[0] << "      "<< sigma_x[0][0] << " "<<sigma_x[0][1] << " "<<sigma_x[0][2] << " "<< sigma_x[0][3] << "                   "<< du[0] << "         "<< lu[0] << "         "<< u2_lu[0]<<"        -------" <<endl;
                cout << "  "<<u[1] << "         "<< sigma_x[1][0] << " "<<sigma_x[1][1] << " "<<sigma_x[1][2] << " "<< sigma_x[1][3] << "                   "<< du[1] << "         "<< lu[1] << "         "<< u2_lu[1]<<"        -------" <<endl;
                cout << "  "<<u[2] << "         "<< sigma_x[2][0] << " "<<sigma_x[2][1] << " "<<sigma_x[2][2] << " "<< sigma_x[2][3] << "                   "<< du[2] << "         "<< lu[2] << "         "<< u2_lu[2]<<"        --------" <<endl;
                cout << "  "<<u[3] << "         "<< sigma_x[3][0] << " "<<sigma_x[3][1] << " "<<sigma_x[3][2] << " "<< sigma_x[3][3] << "                   "<< "--------- don't care -----------------" <<endl;
                cout << "-----------------------------------------------------------------------------"<< endl;
                
                #pragma endregion cout_for_checking_Berlekamps
                
            }

            #pragma region error_count
            for(int i = 0; i < info_len; i++)
            {
                decode_polynomial[i] = y_hat[ i + coding_len - info_len]; //因為 systemaric form, information 在最後k個bits
            }

            int count = 0;
            for (int i = 0; i < info_len; i++)
            {
                if(decode_polynomial[i] != message[i])
                {
                    count = count + 1; 	
                }
            }
            error = error + count;
            #pragma endregion error_count

            if((SNR_dB_2 > 10) && (times > 3000000))
                break;

        }
        //cout << (1.0*error)/(times*info_len) << endl;
        cout << "SNR : " <<SNR_dB_2 << "   BER : " << (1.0*error)/(times*info_len) << endl;
    }
    return 0;
}


int find_symdrome(vector<int> y_hat, int alpha_power, vector<int> G) //find symdrome with decimal term, ex：1+a+a^3 = (a3 a2 a1 a0) = (1 0 1 1) = 11
{
    int symdrome = 0;
    vector<int> alpha_temp;

    for (int i = 0; i < coding_len; i++)
    {
        if(y_hat[i] == 1)
            alpha_temp.push_back( (i * alpha_power) % coding_len );

        /*
        for (int i = 0; i < alpha_temp.size(); i++) 
            cout << alpha_temp[i] <<" ";
        cout << endl;
        */
    }

    for (int i = 0; i < alpha_temp.size(); i++) 
    {
        //cout << symdrome <<"^"<< GF_table_decimal[alpha_temp[i]]<<"=";
        symdrome = symdrome ^ GF_table_decimal[alpha_temp[i]];
        //cout << symdrome <<" "<<endl;
    }
    
    return symdrome;
}

int find_next_du(int current_step, vector<int> sigma_x, vector<int> symdrome_power, vector<int> GF_table_decimal, int lu)
{
    vector<int> dtmp;
    int du_next_decimal = 0;
    int du_next_power;

    //cout << "du+1 = ";
    for (int q = 0; q <= lu; q++)
    {
        //cout << "s:  " << symdrome_power[2*(current_step-1)+3-q-1] <<"   sigma: " << sigma_x[q] <<endl;
        if(symdrome_power[2*(current_step-1)+3-q-1] == 15 || sigma_x[q] == 15 )
            dtmp.push_back( 15 );
        else
            dtmp.push_back( (symdrome_power[2*(current_step-1)+3-q-1] + sigma_x[q]) % 15 );
        //cout << dtmp[q] <<" " << endl;
    }
    //cout << endl;

    
    for (int i = 0; i < dtmp.size(); i++)
    {
        //cout << du_next_decimal << " ^ "<< GF_table_decimal[dtmp[i]] <<" = ";
        du_next_decimal ^= GF_table_decimal[dtmp[i]];
        //cout << "du_next_decimal : "<< du_next_decimal <<" " << endl ;
    }

    for (int j = 0; j <= coding_len; j++)
    {
        if (du_next_decimal == GF_table_decimal[j])
        {
            du_next_power = j;
        }
    }
    //cout << du_next_power << endl ;
    
    return du_next_power;
}

int find_rho(int current_step, vector<int> u2_lu, vector<double> u, vector<int> du)
{
    int rho = -99;
    for (int i = (current_step-1); i >= 0; i--)
    {
        if( du[i] == 15 )
        {
            continue;
        }
        else
        {
            if(u2_lu[i] > rho)
                rho = i;
            else
                rho = rho;
        }
    }
    return rho;
}

vector<int> find_next_sigma_x(int current_step, int rho, vector<vector<int>> sigma_x, vector<double> u, vector<int> du,
                                vector<int> symdrome_power, vector<int> GF_table_decimal)
{
    int X_shifttimes, X_scale_power;
    vector<int> sigma_x_next(minimun_poly_order,0);
    vector<int> sigma_x_rho_shift(minimun_poly_order,0);
    vector<int> sigma_x_next_decimal(minimun_poly_order,0);

    for (int i = 0; i < sigma_x_rho_shift.size(); i++)
    {
        sigma_x_rho_shift[i] = sigma_x[rho][i];
    }
        

    X_shifttimes = int( 2*( (current_step-1.0) - u[rho] ) );
    X_scale_power = (du[current_step] - du[rho]) % 15;
    if(X_scale_power < 0)
        X_scale_power += 15;
    sigma_x_rho_shift.insert(sigma_x_rho_shift.begin(), X_shifttimes, 15);
    for (int i = 0; i < X_shifttimes; i++)
    {
        sigma_x_rho_shift.pop_back();
    }
    /*
    for (int i = 0; i < sigma_x_rho_shift.size(); i++)
    {
        cout <<  sigma_x_rho_shift[i] << " ";
    }
    cout << endl;
    */

    for (int i = 0; i < sigma_x_rho_shift.size(); i++)
    {
        if (sigma_x_rho_shift[i] == 15)
        {

        }
        else
        {
            sigma_x_rho_shift[i] += X_scale_power;
            sigma_x_rho_shift[i] = sigma_x_rho_shift[i] % 15;
        }
        //cout <<  sigma_x_rho_shift[i] << " ";
    }
    //cout << endl;


    for (int i = 0; i < sigma_x_rho_shift.size(); i++)
    {
        sigma_x_next_decimal[i] = GF_table_decimal[sigma_x[current_step][i]] ^ GF_table_decimal[sigma_x_rho_shift[i]];
        /*
        if (sigma_x_next_decimal[i] > 1000 || sigma_x_next_decimal[i] < -1000)
        {
            cout << GF_table_decimal[sigma_x[current_step][i]] << " " <<endl;
            cout << GF_table_decimal[sigma_x_rho_shift[i]] << " "<<endl;
        }
        */

        //cout <<  sigma_x_next_decimal[i] << " ";
    }
    //cout << endl;


    for (int i = 0; i < sigma_x_next_decimal.size(); i++)
    {
        for (int j = 0; j <= coding_len; j++)
        {
            if (sigma_x_next_decimal[i] == GF_table_decimal[j])
            {
                sigma_x_next[i] = j;
            }
        }
        //cout <<  sigma_x_next[i] << " ";
    }
    //cout << endl;


    
    return sigma_x_next;


}
