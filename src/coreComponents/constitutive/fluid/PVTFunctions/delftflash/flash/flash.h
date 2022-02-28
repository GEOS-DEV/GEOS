#pragma once

class Flash {
protected:
	double p, T;
    int NC, NP;
	std::vector<int> ci, pi, eos;
	std::vector<double> z, V, x;

public:
	Flash() {};
	virtual bool runFlash(double pres, double Temp, std::vector<double> zc, std::vector<std::string> comp, std::vector<std::string> ph, std::vector<int> phase_eos) = 0;
	
	std::vector<double> &getV() { return V; }
	std::vector<double> &getx() { return x; }
	std::vector<std::string> getPhases();
	
protected:
	void getIndices(std::vector<std::string> components, std::vector<std::string> phases);
	void normalize0();

};

class SSIFlash : public Flash {
public:
	bool runFlash(double pres, double Temp, std::vector<double> zc, std::vector<std::string> comp, std::vector<std::string> ph, std::vector<int> phase_eos) override;
};

class SSIHydrateFlash :public Flash
{
private:
	std::vector<double> f, x_H, f_wH;
	double f_w;
	std::vector<bool> H_phases;
	int water_index;
	int H_sum{ 0 }; // counts how many hydrate phases are stable
	bool equilibrium = false;
	
public:
	bool runFlash(double pres, double Temp, std::vector<double> zc, std::vector<std::string> comp, std::vector<std::string> ph, std::vector<int> phase_eos) override;
	// bool runEquilibriumFlash();

	// Kinetic conditions output
	double &getfw() { return f[water_index]; }
	std::vector<double> &getfwH() { return f_wH; }
	std::vector<double> &getxH() { return x_H; }

};
