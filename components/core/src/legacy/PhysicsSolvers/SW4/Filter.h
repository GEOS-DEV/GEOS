#ifndef SW4_FILTER_H
#define SW4_FILTER_H

#include <vector>
#include <iostream>

   enum FilterType{lowPass, bandPass};

   class Filter
   {
     friend std::ostream& operator<<(std::ostream& output, const Filter& s);
   public:
      Filter(FilterType type, unsigned int numberOfPoles, unsigned int numberOfPasses, double f1, double f2);
      ~Filter();
      void evaluate(int N, double *u, double *mf);
      void computeSOS(double dt);
      double estimatePrecursor();
      //      std::ostream& operator<<( std::ostream& output );
      FilterType get_type(){return m_type;}
      int get_order(){return m_poles;}
      int get_passes(){return m_passes;}
      double get_corner_freq1(){return m_f1;}
      double get_corner_freq2(){return m_f2;}

      class Polynomial
      {
      	 friend std::ostream& operator<<(std::ostream& output, const Polynomial& s);
      public:
	 Polynomial(){m_c[0]=m_c[1]=m_c[2]=0;}
	 Polynomial( double c[3] ){m_c[0]=c[0];m_c[1]=c[1];m_c[2]=c[2];}
	 double coeff( unsigned int q ){return m_c[q];}
	 double m_c[3];
      };
      class SecondOrderSection
      {
      public:
	 SecondOrderSection(Polynomial &num, Polynomial &den )
	 {  for( int q=0 ; q < 3 ; q++) {m_n.m_c[q]=num.m_c[q];m_d.m_c[q]=den.m_c[q];} }
	 double numer(unsigned int q) const { return m_n.m_c[q];}
	 double denom(unsigned int q) const { return m_d.m_c[q];}
	 Polynomial m_n, m_d;
      private:   
	 SecondOrderSection(){m_n.m_c[0]=m_n.m_c[1]=m_n.m_c[2]=m_d.m_c[0]=m_d.m_c[1]=m_d.m_c[2]=0;};
      };

   private:   
      Filter();
      double realPoleBP(double f1, double f2, double dt, SecondOrderSection *&sos_ptr);
      double complexConjugatedPolesBP(double f1, double f2, double dt, double alpha, 
				      SecondOrderSection *&sos1_ptr, SecondOrderSection *&sos2_ptr);
      double realPoleLP(double fc, double dt, SecondOrderSection *&sos_ptr);
      double complexConjugatedPolesLP(double fc, double dt, double alpha, SecondOrderSection *&sos_ptr);
      void a2d(double n[3], double d[3], Polynomial &b, Polynomial &a);

      FilterType   m_type;
      double       m_dt, m_f1, m_f2;
      unsigned int m_passes, m_poles;
      std::vector<SecondOrderSection*> m_SOSp;
      int    m_numberOfSOS;
      int    m_real_poles, m_complex_pairs;
      bool   m_initialized;
      double m_pole_min_re;
   };

#endif
