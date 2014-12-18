extern "C" {
	typedef struct {
		union {
			double data[2];
			double real, imag;
		};
	} mycomplex;
}
