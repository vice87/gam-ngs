#ifndef _NUCLEOTIDE_
#define _NUCLEOTIDE_

//! Enum representing the type of a base.
typedef enum 
{ 
    A, /*!< Adenine */
    T, /*!< Thymine */
    C, /*!< Cytosine */
    G, /*!< Guanine */
    N, /*!< Any of the four types. */
    LAST_BASE
} __attribute__((packed)) BaseType;

//! Class implementing a nucleotide.
class Nucleotide 
{
private:
    BaseType _base; //!< base type.
  
public:
    Nucleotide();
    Nucleotide(const Nucleotide& orig);
    Nucleotide(const BaseType base_val);
    Nucleotide(const char base);
   
    const BaseType&
    base() const;
 
    const Nucleotide&
    operator=(const Nucleotide& orig); 

    const Nucleotide&
    operator=(const char orig); 

    bool
    operator==(const Nucleotide& a); 

    bool
    operator!=(const Nucleotide& a); 

    operator char() const;

    friend Nucleotide complement(const Nucleotide &orig);
};

Nucleotide complement(const Nucleotide &orig);
std::ostream& operator<<(std::ostream& os, const Nucleotide& base); 

#endif // _NUCLEOTIDE_
