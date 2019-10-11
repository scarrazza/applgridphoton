#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>
#include "fastNLOTools.h"

using namespace std;
using namespace say;
using namespace fastNLO;

namespace fastNLOTools {

   //________________________________________________________________________________________________________________ //
   bool CheckVersion(int version ){
      if ( fastNLO::CompatibleVersions.count(version) == 0 ) {
	 error["fastNLOTools::CheckVersion"]<<"This table version ("<<version<<") is incompatible with this fastNLO code."<<endl;
	 error["fastNLOTools::CheckVersion"]<<"Supported table versions are:";
	 for ( auto i : fastNLO::CompatibleVersions ) error>>" "<<i;
	 error>>""<<endl;
	 error["fastNLOTools::CheckVersion"]<<"Exiting."<<endl;
	 exit(1);
	 return false;
      }
      return true;
   }


   //________________________________________________________________________________________________________________ //
   int ReadVector(vector<double >& v, istream& table , double nevts ){
      //! Read values according to the size() of the given vector
      //! from table (v2.0 format).
      for( unsigned int i=0 ; i<v.size() ; i++){
	 table >> v[i];
	 v[i] *= nevts;
	 if ( !isfinite(v[i]) ) {
            error["ReadVector"]<<"Non-finite number read from table, aborted! value = " << v[i] << endl;
            error["ReadVector"]<<"Please check the table content." << endl;
            exit(1);
         }
      }
      return v.size();
   }

   //________________________________________________________________________________________________________________ //
   int ReadUnused(istream& table ){
      //! Read values, which are not known to the current code.
      int nLines = 0;
      table >> nLines;
      if ( nLines==fastNLO::tablemagicno ) {
	 error["ReadUnused"]<<"Number of lines identical to magic number. Exiting."<<endl; exit(3);
      }
      string sUnused;
      if ( nLines > 0 ) std::getline(table,sUnused); // discard empty space due to precendent >>
      for( int i=0 ; i<nLines ; i++)
	 std::getline(table,sUnused);
      return nLines;
   }


   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<double >& v, istream& table , int nProcLast , double nevts ){
      int nn = 0;
      if ( nProcLast == 0 ) {
         table >> nProcLast; 
	 nn++;
      }
      v.resize(nProcLast);
      for(unsigned int i0=0;i0<v.size();i0++){
	 table >> v[i0];
	 v[i0] *= nevts;
	 nn++;
	 if ( !isfinite(v[i0]) ) {
            error["ReadFlexibleVector"]<<"Non-finite number read from table, aborted! value = " << v[i0] << endl;
            error["ReadFlexibleVector"]<<"Please check the table content." << endl;
            exit(1);
         }
      }
      return nn;
   }

   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<unsigned long long >& v, istream& table , int nProcLast , double nevts ){
      int nn = 0;
      if ( nProcLast == 0 ) {
         table >> nProcLast; 
	 nn++;
      }
      v.resize(nProcLast);
      for(unsigned int i0=0;i0<v.size();i0++){
         char buffer[256];
         table >> buffer;
         double value = atof(buffer);
	 v[i0]  = value;
	 v[i0] *= nevts;
	 nn++;
      }
      return nn;
   }


   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<std::string >& v, istream& table , int size , double nevts ){
      if ( size == 0 ) table >> size;
      v.resize(size);
      if ( size > 0 ) std::getline(table,v[0]); // discard empty space due to precendent >>
      for( auto& i : v) {
	 std::getline(table,i);
      }
      return v.size() + 1;
   }


   //________________________________________________________________________________________________________________ //
   int ReadFlexibleVector(vector<int >& v, istream& table , int size , double nevts ){
      if ( size == 0 ) table >> size;
      v.resize(size);
      for( auto& i : v) {
	 table >> i;
      }
      return v.size() + 1;
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v7d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3, dim4, dim5, dim6 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v6d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3, dim4, dim5 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }




   //________________________________________________________________________________________________________________ //
   void ResizeVector( v5d& v, int dim0 , int dim1, int dim2, int dim3, int dim4 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3, dim4 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v4d& v, int dim0 , int dim1, int dim2, int dim3 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2, dim3 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v3d& v, int dim0 , int dim1, int dim2 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1, dim2 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v2d&  v, int dim0 , int dim1 ){
      if ( dim0 > 0 ){
         v.resize(dim0);
         for ( int i= 0 ; i<dim0 ; i++)
            ResizeVector( v[i] , dim1 );
      } else {
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }


   //________________________________________________________________________________________________________________ //
   void ResizeVector( v1d& v, int dim0 ){
      if ( dim0 > 0 )
         v.resize(dim0);
      else{
         error["fastNLOTools::ResizeVector"]<<"Cannot resize table, because dimension is <= zero (dim0="<<dim0<<"). Exiting"<<endl;
         exit(1);
      }
   }



   //______________________________________________________________________________
   void ResizeFlexibleVector(vector<double >& v, const vector<double >& nom) {
      v.resize(nom.size());
   }

   //________________________________________________________________________________________________________________ //
   void PutBackMagicNo(istream& table){
   //! Put magic number back
      for(int i=0;i<(int)(log10((double)tablemagicno)+1);i++){
         table.unget();
      }
      table.unget();
   }

   //______________________________________________________________________________
   bool ReadMagicNo(istream& table) {
      //! read and crosscheck magic number
      if (table.eof()){
	 error["ReadMagicNo"]<<"Cannot read from file. Exiting"<<endl; 
	 exit(3);
      }
      string line;
      std::getline(table,line);
      if ( line=="" ) std::getline(table,line);  // last one was '<<'
      if( line != std::to_string(tablemagicno)){
         error["ReadMagicNo"]<<"Found '"<<line<<"' instead of "<<tablemagicno<<"."<<endl;
	 error["ReadMagicNo"]<<"Did not find magic number, aborting!"<<endl;
	 error["ReadMagicNo"]<<"Please check compatibility of tables and program version. Exiting."<<endl;
	 exit(2);
         return false;
      };
      return true;
   }


} // end namespace fastNLO
