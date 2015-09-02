//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#ifndef SIZED_ARRAY_H
#define SIZED_ARRAY_H
#include <algorithm>
/** This class manages arrays as if they were STL vectors. This is useful to interact with Geant4 methods which require arrays instead of vectors.
 */

template <typename V> class sized_array {
  public:
    ///Constructor of an array of size "size"	 
    sized_array (int newsize) { 
        _size(newsize);
	_ptr = new V[_size]; 
    }

    ///Initialization of an empty array
    sized_array (){
    _size=0;
    _ptr=NULL;
    }

    ///clear of an array
    void clear(){
	    _size=0;
	    _ptr=NULL;
    }

    ///Destructor
    ~sized_array () { delete [] _ptr; }
    
    ///Copy constructor
    sized_array (const sized_array<V>& r): _size(r._size) { 
	    _ptr = new V[_size]; 
	    std::copy (r._ptr, r._ptr + _size, _ptr); }
   
   ///Vector allocation if it is still empty
    void allocate (int newsize) {
	   if(!_ptr) {
		_size = newsize ; 
		_ptr = new V[_size];}
    } 
    
    ///Vector reallocation according to the new size "size"
    void reallocate (int newsize) { 
	V* new_ptr = new V[newsize]; 
	if(_ptr) { 
	std::copy (_ptr, _ptr + std::min(_size, newsize), new_ptr);
	} 
	_size = newsize; 
	_ptr = new_ptr; 
    }

    ///Redefinition of the assignment operator
    const sized_array<V>& operator= (const sized_array<V>& r)  { 
	  delete [] _ptr; 
	 _size = r._size;  
	 _ptr = new V[_size]; 
	 std::copy (r._ptr, r._ptr + _size, _ptr); 
    }

    ///It returns the pointer of the object
    operator V* () { return _ptr; } 
    
    operator const V* () const { return _ptr; } 
    
    ///Access to the i-th element
    V& operator[] (int i) { return _ptr[i]; }

    ///Alternative way to acces the i-th element
    V at (int i) {return _ptr[i]; }

    const V& operator[] (int i) const { return _ptr[i]; }
    
    int size () const { return _size; }
  
  private:

    V* _ptr;
    int _size;

};

#endif
