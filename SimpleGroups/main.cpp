//
//  main.cpp
//  SimpleGroups
//
//  Created by XuXo on 1/01/14. 
//  Copyright (c) 2014 XuXo. All rights reserved.
//

// simpleGroups.cpp : Defines the entry point for the console application.
// This is a quick and dirty solution for verifying a group of a certain order is simple (ie no non trivial normal subgroups) as I've done way too many of these problems.

//#include "stdafx.h"		//vs2010 dependent
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>

#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <math.h>
#include <map>

#define VERBOSE
#define DEBUG

using namespace std;


/******************************************************************************/
static vector<int> divisors (int num);
static int factorial (int num);
static int reverseFactorial (int num);
static bool trick1 (int num) ;
static bool trick2 (int num);
static bool classEquation (int order, vector<int> &factors);
static bool simpleSylow (int num);
/******************************************************************************/


int main(int argc, const char *argv[])
{
    
	//vector<int>sdfsd = divisors(105);
    //testing some cases
    //trick1(392);        //should return 7;
    //trick1(1452);       //should return 11;
	//trick2(480);		  //should return 5;
    simpleSylow (28);
    
    //file.close();
    cin.get();
    return 0;
}
//720, 756, and 840 would fail miserably

//this name suggests a complexSylow which does exist, a generalized version but it's difficult to implement from a counting perspective
static bool simpleSylow (int num) {
    
    
    std::map<int,int> factors;
    int order_of_group = num;
    cout<<"The prime decomposition of the number is :";
    for (int i = 2; i< order_of_group/2; i++)
    {
        int exp = 0;
        bool is_factor = false;
        
        while (num % i  == 0)
        {
            is_factor = true;
            exp++;
            cout<<i<<" ";
            num = num/i;
        }
        
        if(is_factor)
        {
            map<int, int>::iterator it = factors.find(i);
            if( it == factors.end())
                factors.insert(std::make_pair(i, exp));
        }
    }
	cout<<endl<<endl;
    
    /*-----------------------------test for class equation for centralizers of distinct conjugacy classes-----------------------------*/
    vector<int> list_of_divisors = divisors(order_of_group);
	if(!classEquation(order_of_group, list_of_divisors))
	{
        cout<<"G is non-simple due to the Class Equation. ";
        #ifdef VERBOSE
        cout<<"Remember the Class Equation states that |G| = |Z(G)| + Sum over i of (|G:C_G(gi)|) where g, g2,...gr the representatives of the distinct conjugacy classes of G not contained in the center Z(G). "<<endl;
        cout<<"Assuming G is simple implies Z(G) = 1 so our equation reduces to "<<order_of_group<<" = 1 + Sum over i of (|G:C_G(gi)|)"; ;
        cout<<"Picking any arbitrary divisor of "<<order_of_group<<", we know the RHS has at least one term not divisible by this particular divisor which means it must divide into one of the remaining divisors"<<endl;
        cout<<"ie |G:C_G(a)| = d for some element a in G and divisor d of |G|.  This is equivalent to showing |G||d! since order of simple non-abelian group divides factorial of index of proper subgroup";
        cout<<"(and we know only abelian simple groups are those of order p).  So as long as we find a divisor for which this does not hold, it can't be simple."<<endl<<endl;
        #endif
		return false;
	}
    
    
	/*-----------------------------test for order = 2k where k is odd means non-simple-----------------------------*/
	if(order_of_group % 2 == 0)
    {
		int half = order_of_group/2;
		if(half % 2 == 1)
        {
			cout<<"The order is divisible by 2 but not 4 and therefore has a normal subgroup of index 2 forcing G to be non-simple ";
            #ifdef VERBOSE
            cout<<"due to the following observation."<<endl;
			cout<<"Since all elements of G act on G by translation, we have a homomorphism phi : G -> S_2k.  Denote its image by H.  It's enough to show H has a subgroup of index 2."<<endl;
			cout<<"Since H contains an odd permutation, H is not in A_2k.  From HA_2k = S_2k, we derive that :"<<endl;
			cout<<"|S_2k|/|A_2k| = |H| / |H intersect A_2k|."<<endl;
			cout<<"Hence |H|/|H intersect A_2k| = 2."<<endl<<endl;
            #endif
            return false;
		}
	}
    
    
	/*-----------------------------test for order = p^n for prime p, p groups are always non-simple-----------------------------*/
    if(factors.size() == 1)
    {
        
        //Sylow will force np to be 1 but there's a much easier argumen  to short circuit this calculation
		cout<<"This is exactly a p-group which is known to be non-simple ";
        
        #ifdef VERBOSE
        cout<<"due to the following observation."<<endl;
        cout<<" For any finite group G it is obviously true that for any proper subgroup H, the normalizer N_G(H) is strictly larger than H itself.  Hence for a subgroup H of order p^(n-1) which is sure to exist the normalizer must be equal to all of G, implying H is normal.  In fact a p-group of order p^n has a normal subgroup of order pK for each 0<=k<=n"<<endl;
        #endif
		return false;
	}
    
	
	/*-----------------------------test for order = pqr, a product of 3 primes-----------------------------*/
	//if(factors.size() == 3)
	bool pqr = true;
	if(factors.size() == 3)
    {
		map<int, int>::iterator it;
		for (it = factors.begin(); it != factors.end(); ++it)
			if (it->second != 1)
				pqr = false;
        
		if(pqr)
        {
			cout<<"The order is of the form pqr and is non-simple ";
            
            #ifdef VERBOSE
            cout<<"due to the following observation."<<endl;
			cout<<"Assume WLOG p < q < r.  We will show that G has a normal Sylow subgroup for either p, q, or r.  Applying Sylow's criterion forces the follow."<<endl<<endl;
			cout<<"np in { 1, q, r, qr }"<<endl;
			cout<<"nq in { 1, p, r, pr }"<<endl;
			cout<<"nr in { 1, p, q, pq }"<<endl<<endl;
			cout<<"Moreover, since p < q < r, we have p not congruent to 1 mod q, p not congruent to 1 mod r, q not congruent to 1 mod r.  So narrowing the list further: "<<endl<<endl;
			cout<<"np in { 1, q, r, qr }"<<endl;
			cout<<"nq in { 1, r, pr }"<<endl;
			cout<<"nr in { 1, pq }"<<endl<<endl;
			cout<<"Now suppose G is simple so none of np, nq, nr is equal to 1. This forces nr = pq.  The Syl_r(G) in turn are cyclic and intersect trivially, so G has pq(r-1) elements of order r.  Similarly, G has at least q(p-1) elements of order p and r(q-1) elements of order q. Then the number of elements in G is at least pq(r-1) + q(p-1) + r(q-1) = pqr + rq - r - q. However, because 3 <= q < r, q+r < qr, so that G has too many elements, a contradiction."<<endl;
            #endif
            return false;
		}
    }
    
    
    /*-----------------------------test for order = 4p^n for prime p > 2, realize this is just a special case of Burnside-----------------------------*/
    if(factors.size() == 2)
    {
        if (factors[2] == 2)
        {
            cout<<"The order is of the form 4p^n and is non-simple ";
            
            #ifdef VERBOSE
            cout<<"due to the following observation."<<endl;
            cout<<"First consider when p>= 5. Then p>4 so np must be 1 due to congruence relation. Syl_p(G) must be normal and hence G non-simple"<<endl;
            cout<<"now let p be 3.  Consider n = 1.  Then |G| = 12 and if Syl_p(G) is normal, we are done.Otherwise, G is isomorphic to A4 by a corollary in Dummit & Foote which is known to be non-simple."<<endl;
            cout<<"we deduce n >= 2.  Quoting the theorem that if H is a subgroup of finite index n in a gorup G, then there exists a homomorphism phi: G -> Sn with kernerl in H"<<endl<<endl;
            
            cout<<"We must have a homomorphism phi: S -> S4 with kernel in Syl_p(G).  Since n >= 2, we have |G| > 9*4 = 36 > |S4| = 24.  Therefore phi cannot be injective which means kernel is nontrivial and so G is non-simple"<<endl<<endl;
            #endif
            return false;
        }
    }
    
    
    /*-----------------------------test for order = p^n(p+1) for some prime p-----------------------------*/
    if(trick1(order_of_group))
    {
		cout<<"The order is of the form p^k(p+1) and is non-simple ";
        
        #ifdef VERBOSE
        cout<<"due to the following observation"<<endl;
        cout<<"We suppose there is such a simple group.  Then np is forced to equal p+1 due to Sylow.  Consider the action of G on Syl_p(G) via usual conjugation which gives a homomophism from G -> Sp+1"<<endl;
        cout<<"Since G is simple then either kernel of the map is trivial or all of G.  If kernel was G, then the image of G in the symmetric group is just the identity, so G does not move any of its Syl_(G) around"<<endl;
        cout<<"But we know G acts transitively on the conjugates which is only possible if p+1 = 1, contradiction"<<endl;
        
        cout<<"If Ker(G) = 1, then |G| divies (p+1)! which is possible only when p = 1, contradiction.  Indeed, taking p=2, non abelian group of order 6 has order p(p+1) and has exactly p+1 Syl_p(G) subgroups."<<endl<<endl;
        #endif
        return false;
	}
    
    
    /*-----------------------------test for order = (p^2-p)( p^-1)-----------------------------*/
	if(trick2(order_of_group))
    {
		cout<<"The order is of the form (p^2 - p)(p^2 - 1) which is the order of Gl_2(Fp) with p+1 Sylow P subgroups and thus non-simple ";
        
        #ifdef VERBOSE
        cout<<"due to the following observation"<<endl;
		cout<<"Sylow tells us the number of sylow p-subgroups divides (p-1)^2(p+1) and is congruent to 1 mod p.  It turns out the structly upper triangular matrixes form a sylow P-subgroup (easy to verify) and a simple counting argument tells us there are p(p-1)^2 of these since diagonal entries are arbitrary nonzero field elements while upper right entry is unrestricted."<<endl;
		cout<<"[G: UT_2(Fp)] = [G: N_G(SUT_2(Fp))][N_G(SUT_2(Fp)): UT_2(Fp)] where SUT denotes strictly upper triangular. "<<endl;
		cout<<"By Sylow's Theorem, then, np divies p+1.  But since np is congruent to 1 mod p, either np = 1 or np = 1+p.  np is clearly not 1 since strictly upper and strictly lower triangular matrices are distinct syloy p-subgroups of GL_2(Fp)"<<endl;
        #endif
		return false;
	}
    
    
    /*-----------------------------test for Burnside very straightforward-----------------------------*/
    if(factors.size() == 2)
    {
		cout<<"Order of the group is of the form p^a*q^b and therefore due to Burnside's and thus non-simple ";
        
        #ifdef VERBOSE
        cout<<"due to the following observation.  This is because all groups of this form must be solvable and only finite solvable simple groups are the cyclic groups of prime order."<<endl;
        
        #endif
		return false;
	}
    
    
    /*-----------------------------test for order = p1 p2 p3... where pi's are distinct and pi does not divide pj-1 for all i,j then G is nonsimple, infact it's cyclic (this is too involved though)-----------------------------*/
	vector<int>primes;				//work off a vector instead of a map since we need fast access on multiple iterative searches
	bool power_free = true;			//(in the spirit of square free)
	map<int, int>::iterator it;
	for (it = factors.begin(); it != factors.end(); ++it)
    {
		if( it->second == 1)
			primes.push_back(it->first);
		else power_free = false;
    }
    
	bool divisibility = false;
	if(power_free)
    {
		for (int i = 0; i < primes.size(); i++)
			for (int j = i+1; j < primes.size(); j++)
				if ((primes[j] % (primes[i]-1) != 0) || (primes[i] % (primes[j]-1) != 0))		// p_i does not divide p_j and p_j does not divide p_i
					divisibility = true;
        
		if(divisibility == false)
        {
			cout<<"The order is of the form p1 p2 p3...pr  where pi's are distinct and p_i does not divide p_j-1 for all i,j G is therefore non-simple ";
            
            #ifdef VERBOSE
            cout<<"due to the following observation"<<endl;
			cout<<"We proceed by induction on the width of G. (Recall that the width of a finite group is the number of prime divisors (including multiplicity) in its prime factorization)."<<endl;
			cout<<"For the base case, suppose G has width 1. Then |G| = p, and clearly p does not divide p-1. Moreover, every group of order p (a prime) is cyclic. "<<endl;
			cout<<"For the inductive step, suppose that every finite group of width at most k and whose orderÌs prime divisors are distinct and satisfy the divisibility requirement is cyclic. Let G be a finite group of width k+1, say |G| = p_1 p_2... p_{k+1}, whose orderÌs prime divisors are distinct and satisfy the divisibility requirement."<<endl;
			cout<<"Note that every proper subgroup of G has width at most k and satisfies the distinctness and divisibility criteria; by the induction hypothesis, every proper subgroup of G is cyclic. Using this previous exercise, G is not simple (whether or not it is abelian"<<endl<<endl;
            #endif
			return false;
		}
	}
    
    
	/*-----------------------------test for Feit-Thompson very straightforward-----------------------------*/
	if((factors.size() > 1) && (order_of_group % 2 == 1))
    {
		cout<<"The order of the group is of the form p1^a1 * p2^a2 * ... pj^aj where j>1 and hence non-simple ";
        
        #ifdef VERBOSE
        cout<<"since the only simple groups of odd order are those of prime order.  This follows directly from Feit-Thompson Theorem that every group of odd order is solvable"<<endl;
        #endif
		return false;
	}
    
    
	/*-----------------------------test for alternating group An where n>= 5, very straightforward-----------------------------*/
	int alternating = order_of_group*2;
	if((order_of_group >= 60) && (reverseFactorial(alternating) >= 5) )
    {
		cout<<"This is A"<<reverseFactorial(alternating)<<" and hence simple"<<endl;
        return false;
    }
    
    
    #ifdef DEBUG
    for (it = factors.begin(); it != factors.end(); ++it) {
        cout<< "prime factor : "<<it->first<<" , power :  "<<it->second<<endl;
    }
    #endif
    
	//assuming it passed the basic tests, we will use some more elaborate counting techniques now
    std::map<int,vector<int>> np;
    for (it = factors.begin(); it != factors.end(); ++it)
    {
        
        //cout<<it->first<<" ";
        vector<int> np_list;
        map<int, int>::iterator it_inner;
        for( it_inner = factors.begin(); it_inner != factors.end(); it_inner++)
        {
            //cout<<itInner->first<<" ";
            if (it->first != it_inner->first)
            {
                np_list.push_back(pow((float)it_inner->first, it_inner->second));
                
                //cout<<itInner->first<<" ";
                //cout<<pow((float)itInner->first, itInner->second)<<" ";
            }
        }
        np.insert(std::make_pair(pow((float)it->first, it->second), np_list));
    }
    cout<<endl;
    
    
    #ifdef DEBUG
    map<int, vector<int>>::iterator it2;
    for (it2 = np.begin(); it2 != np.end(); ++it2)
    {
        cout << "n"<<it2->first <<" possible values : ";
        for (int i = 0; i<it2->second.size(); i++)
            cout<<it2->second[i]<<" ";
        cout<<endl;
    }
	cout<<endl;
    #endif
    
    
    std::map<int, vector<int>> NP;
    for (it2 = np.begin(); it2 != np.end(); ++it2)
    {
        vector<int> NPList;
        int product=1;
        for (int i = 0; i<it2->second.size(); i++) {
            product = product*it2->second[i];
        }
        
        vector<int> div = divisors(product);
        
        //a number is technically a divisor of itself, should append this to the divisors() function
        //div.push_back(product);
		cout<<"divisors of m : ";
        for(int i =0; i<div.size(); i++)
            cout<<div[i]<<" ";
        cout<<endl;
        
        for(int i = 0;i< div.size(); i++)
        {
            vector<int> facs2 = divisors(it2->first);
            
            //the first element is always of the form p^a for some a, so divisors should return 1, p, p^2, p^3, ...p^(a-1), p^a always so we need the 2nd element for the sylow modulus criteria
            if( div[i] % facs2[1] == 1)
                
                //cout<<div[i]<<" ";
                NPList.push_back(div[i]);
        }
        NP.insert(std::make_pair(it2->first, NPList));
    }
    
    #ifdef DEBUG
    cout<<endl<<endl;
    map<int, vector<int>>::iterator it3;
    for (it3 = NP.begin(); it3 != NP.end(); ++it3)
    {
        cout << "n"<<it3->first <<" filtered list of possible values : ";
        for (int i = 0; i<it3->second.size(); i++)
            cout<<it3->second[i]<<" ";
        cout<<endl;
    }
    cout<<endl;
    #endif
    
    
    /*-----------------------------test for number of sylow groups being 1-----------------------------*/
    
    //if any one of these is is only 1 then we definitely have a normal subgroup hence not simple
    for (it3 = NP.begin(); it3 != NP.end(); ++it3)
    {
        if(it3->second.size() == 1)
        {
            cout<<"n"<<it3->first<<" can only be 1 and due to Sylow's Theorem, Syl_"<<it3->first<<"(G) must be normal and hence non-simple"<<endl;
            return false;
        }
    }
    cout<<endl;
    
	/*-----------------------------test for |G| divides np!-----------------------------*/
    
	//look at G acting on the set of sylow subgroups and derive contradiction
	//ex. |G|=224, working through sylow tells us n2 must be 1 or 7 and n7 must be 1 or 8, look at the Sylow 2 subgroups, of order 32.
	//suppose n2 is 7, ie 7 subgroups in G.  Let G act on this set of Sylow subgroups via conjugation, we get a homomorphism G -> S7 which must be injective if G is simple since
	//otherwise it would have nontrivial kernel and the kernel of the map is always normal subgroup of G.  But G can't possibly embed in S7 since |G| does not divide 7!
    
	//so we can very easily check this condition to filter down the list of n_ps, redo this to actually trim the list
	for (it3 = NP.begin(); it3 != NP.end(); ++it3)
    {
		bool trivial_kernl = true;
		vector<int> possible_np;
		for (int i = 1; i < (it3->second).size(); i++)
        {		//start at index 1 since index 0 = 1 and its trivial
			if (factorial((it3->second)[i]) % order_of_group == 0)
            {
				possible_np.push_back((it3->second)[i]);
				cout<<"order of the group does divide evenly into "<<(it3->second)[i]<<"factorial"<<endl;
			}
		}
		cout<<endl;
		if(possible_np.empty())
			return false;
	}
    
    
    //since all the sylow groups where exponent is 1 are cyclic they must intersect trivially therefore we can calculate the number of elements of their respective orders and try to
	//get a contradiction if they exceed total number of elements
    int total_elem = 0;
    cout<<"number of groups"<<endl;
    for (it3 = NP.begin(); it3 != NP.end(); ++it3)
    {
        vector<int> fac = divisors(it3->first);
        
        //cout<<fac.size();
        if (fac.size() == 2)
        {       //ie not prime, we are only interested in elements where the exponent is 1
            
            total_elem += (fac[1]-1) * it3->second[1];        //we always pick the smallest np index
            
            //Syl_p(G) is of order p and thus isomorphic to Zp so every non-identity element has order p-1, but this argument is false if it's a Syl_p^k(G) group since there is no
			//reason why the non-identity elements would have order p^k as well
            cout<<" There are "<<it3->second[1]<<" groups of order "<<fac[1]-1<<endl;
            
        }
    }
    
    
    
    /*-----------------------------test for total number of elements exceeding |G| by counting argument-----------------------------*/
    if (total_elem > order_of_group)
    {
        cout<<"There are more elements of finite order than the order of the group itself forcing np to be 1 and hence non-simple"<<endl;
        return false;
    }
    
    
	/*-----------------------------test for something similar, illustrated with an ex.-----------------------------*/
    
	//|G|= 132 = 2^2 * 3 * 11, applying Sylow we arrive at the following list of possible indices
	//n2 is (1,3,11 or 33).  n3 in (1,4,22) n11 in (1,12)
	//if we pick 12 for n11 and 4 for n3 we get 128 elements total which leaves exactly 4 elements which must be the unique sylow 2 subgroup, forcing normality and hence nonsimple
    
	//2x2. |G| = 992 = 2^5 * 31 from Sylow we get either 1 or 32 for n31.  Assume simple so n31 = 32, this gives 32*(31-1) = 960 elements. leaving 32 elements left, which is exactly
	//the number of elements in the Syl_2^5(G) subgroup, hence n2 = 1 and is normal.
    
	cout<<"total elem"<<total_elem<<endl;
	for (it3 = NP.begin(); it3 !=NP.end(); ++it3)
    {
		vector<int> fac = divisors(it3->first);
        
        //cout<<fac.size();
        if (fac.size() != 2)                    //it's a syl p^a group where a is not 1
			total_elem += it3->first;
	}
	if (total_elem == order_of_group) {
		cout<<"There are exactly enough elements left that give 1 unique normal sylow subgroup and hence non-simple"<<endl;
		return false;
	}
	
	/*-----------------------------Finally Classfication-----------------------------*/
    
	//classification theorem, eventually will add this but it would be nice to be able to derive these without have to resort to classifying.  This is definitely easier and more definitive
	//since there are an inexhaustible list of tricks and counting argements.  Basic idea behind the classification approach is that We know finite simple groups are classified as being in
	//one of 19 families, or being one of 26 exceptions
	//Zp, p prime, An (alternating) with n >= 5, one of 16 groups of lie type, one of 26 exceptions, the sporadic groups, of which 20 are subgroups or subquotients of the monster group
    
    
	//generalized Sylow
	//Let P1, P2, P3,...,Pk be sylow -subgroups of G.  if [Pi: Pi intersect Pj] >= p^d whenever i not equal j, then np congruent to 1 mod p^d.  ie notice taking d = 1 gives us Sylow's exactly.
	//The contrapositive is useful in obtaining a contradiction, if np is not congruent to 1 mod p^2, then there exist distinct Syl_p(G) P,Q s.t. |P: P intersect Q| < p^2 and |Q: P intersect Q| < p^2
	//which in effect says |P: P intersect Q| = p and |Q: P intersect Q| = p, and thus P intersect Q is normal in both P and Q (since any subgroup of index p in a p-group is normal).
	//Therefore P intersect Q is normal in <P,Q> as well.  Now |<P,Q>| > |Q| of course since P and Q are distinct so if |G| = p^k * q we know |<P,Q>| must equal |G| and hence <P,Q> = G so P
    //intersect Q is normal in G.  Interestingly enough notice if |G| is of this form, we would have conclude non-simple by Burnside already.
    
    cout<<"Inconclusive...hope I never see this but I will"<<endl;
	return true;
}

//finding all divisors, not just the prime ones
static vector<int> divisors (int num) {
    vector<int> divisors;
    for (int i = 1; i <= num/2; i++)
    {
        if(num % i ==0)
            divisors.push_back(i);
    }
    divisors.push_back(num);
    return divisors;
}

static int factorial (int num) {
	int prod = 1;
	while (num > 1)
    {
		prod = prod*num;
		num--;
	}
	return prod;
}

//checking if a number is a factorial
static int reverseFactorial (int num) {
	int counter = 1;
	while(num > 1)
    {
		if(num % counter != 0)
			return 0;
		num = num / counter;
		counter++;
	}
	cout<<"yes it's"<<counter<<" factorial"<<endl;
	return counter;
}


static bool trick1 (int num) {
    /*testing to see if it's of the form p^n(p+1) for some n
     realize this is roughly p^(n+1) so we iterate through all n and use the n+1 th root of num as an initial estimate to further improve upon
     p^(1/n) = exp( ln (p^(1/n))) = exp ( ln(p)/n))
     
     so how high can n go, well the smallest prime is 2 so our upper bound is the smallelst n such that 2^n > num
     */
    
    float max_nthpower = log((float)num)/log((float)2);
    cout << "highest nth root is we need to consider is " << max_nthpower << endl;
    for (int i = 2; i < max_nthpower; i++)
    {
        
		//we need log(float()) because this is a Visual Studio issue otherwise will complain about ambiguity in overloading
        float nth_root = exp(log((float)num) / i);
        cout<< "currently testing for the "<<i<<"th root using approximation "<<nth_root<<endl;
        
        //we use a +/- 2 cushion just to be safe, this might go negative but that's not a problem really
        for( int j = i - 2 ; j < i + 2; j++)
        {
            cout<<"j is "<<j<<endl;
            if((pow((float) floor(nth_root), j) * (floor(nth_root) + 1)) == num)
            {
                cout<<"the p in questions is "<< floor(nth_root)<<endl;
                return true;
            }
        }
    }
    return false;
}


static bool trick2 (int num) {
    /*testing to see if it's of the form (p^2-1)(p^2-p) for some n
     this is roughly (p-1)(p+1)(p-1)p so roughly p^4 and so we test some 4th roots
	 this is a lot easier than trick1 since we know the order of the root at least
     */
    
    float max_nthpower = log((float)num)/log((float)2);
    cout << "highest nth root is we need to consider is " << max_nthpower << endl;
    float fourth_root = exp(log((float)num) / 4);
    cout<< "currently testing for the "<<fourth_root<<"th root"<<endl;
    
    //we use a +/- 5 cushion just to be safe, this might go negative but that's not a problem really
    for( int j = floor(fourth_root) - 5 ; j < floor(fourth_root) + 5; j++)
    {
        cout<<"j is "<<j<<endl;
		if(pow((float)(j-1), 2) * (j+1) * j == num)
        {
			cout<<"the p in question is "<< j<<endl;
			return true;
		}
    }
    return false;
}


static bool classEquation (int order, vector<int> & factors) {
	//we just need to find a single divisor for which the the class equation fails to hold
	for(int i = 1; i< factors.size(); i++) {
		cout<<"|G|/|C_G(a)| must not be divisible by "<<factors[i]<<" for some a in G"<<endl;
		bool divisible = false;
		for (int j = 1; j<factors.size()-1; j++)
		{
			//so we iterate through all the other divisors and see if it's a viable choice for index G/C_G(a)
			if ((j != i) && (factors[j] % factors[i] != 0))
			{
				cout<<"divisor is "<<factors[j];
				cout<<".  factorial is"<<factorial(factors[j])<<endl;
				if(factorial(factors[j]) % order == 0)
				{
					//as soon as we find one then this eliminates the possibility of disproving the class equation we move onto the next divisor
					cout<<"the order of G divides "<<factors[j]<<"!"<<"so we can move on"<<endl<<endl;
					divisible = true;
					break;
				}
			}
		}
		if (!divisible) cout<<endl;
		if (divisible == false) {
            
			cout<<"We know at least one term (|G|/|C_G(a)| for some a in G) the RHS is not divisible by "<<factors[i]<<"which means it must be equal to one of the remaining divisors in factors[]";
			cout<<"or |G| divides factor[i]! for some i but that is not the case.  In more detail..."<<endl<<endl;
			return false;
		}
	}
}
