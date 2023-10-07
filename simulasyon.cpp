
#include <iostream> 
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <queue>
#include <thread>
#include <chrono>
class Environment {
public:
    Environment(double initialPollutionLevel) : pollutionLevel(initialPollutionLevel) {}

    void emitPollution(double emissionRate) {
        pollutionLevel += emissionRate;
    }

    void simulate(double time, double emissionRate) {
        for (double t = 0; t < time; t += 1.0) {
            emitPollution(emissionRate);
            std::cout << "Zaman: " << t << " saat, Kirlilik Seviyesi: " << pollutionLevel << " birim" << std::endl;
        }
    }
private:
    double pollutionLevel; 
};
class SurfaceTensionSimulator {
public:
    SurfaceTensionSimulator(double surfaceTensionCoefficient)
        : surfaceTensionCoefficient(surfaceTensionCoefficient) {}
        
    double calculateSurfaceTension(double radius) {
        double pi = acos(-1.0);
        return 2.0 * pi * radius * surfaceTensionCoefficient;
    }

private:
    double surfaceTensionCoefficient; 
};
class HydraulicSystem {
public:
    HydraulicSystem(double pistonArea, double loadForce, double initialPosition)
        : pistonArea(pistonArea), loadForce(loadForce), position(initialPosition) {}

    void simulate(double timeStep, double simulationTime) {
        double currentTime = 0.0;
        while (currentTime < simulationTime) {
            double pressure = loadForce / pistonArea;
            double velocity = pressure / pistonArea;  
            position += velocity * timeStep;
            std::cout << "Zaman: " << currentTime << " s, Pozisyon: " << position << " m" << std::endl;
            currentTime += timeStep;
        }
    }
private:
    double pistonArea;   
    double loadForce;    
    double position;     
};
class DCMotorSimulation {
public:
    DCMotorSimulation(double Kp, double Ki, double Kd)
        : Kp(Kp), Ki(Ki), Kd(Kd), targetSpeed(0), currentPosition(0), integral(0), lastError(0) {}

    void setTargetSpeed(double speed) {
        targetSpeed = speed;
    }

    void update(double dt) {
        double error = targetSpeed - currentPosition;
        integral += error * dt;
        double derivative = (error - lastError) / dt;

        double controlSignal = Kp * error + Ki * integral + Kd * derivative;
        currentPosition += controlSignal * dt;

        lastError = error;
    }

    double getCurrentPosition() const {
        return currentPosition;
    }

private:
    double Kp, Ki, Kd;
    double targetSpeed;
    double currentPosition;
    double integral;
    double lastError;
};
class Vehicle {
public:
    double position;
    double velocity;
    double acceleration;

    Vehicle() : position(0.0), velocity(0.0), acceleration(0.0) {}

    void update(double timeStep) {
        velocity += acceleration * timeStep;
        position += velocity * timeStep;
    }
};
class Material {
public:
    double thermalConductivity; 
    double specificHeat;        

    Material(double conductivity, double heat) : thermalConductivity(conductivity), specificHeat(heat) {}
};
class Part {
public:
    double initialTemperature;
    double length;

    Part(double temp, double len) : initialTemperature(temp), length(len) {}
};
double calculateHeatTransfer(double conductivity, double specificHeat, double initialTemp, double finalTemp, double length) {
    double deltaTemp = finalTemp - initialTemp;
    return (conductivity * deltaTemp) / length * specificHeat;
}

class Akiskan {
public:
    double Hiz;       
    double Yogunluk;  
};
double HesaplaBasinc(const Akiskan& akiskan) {
    double sabit = 1.0;
    double basinc = sabit - 0.5 * akiskan.Yogunluk * akiskan.Hiz * akiskan.Hiz;
    return basinc;
}
int main()
{
    int cho;
    std::cin >> cho;
    if (cho == 1)
    {
        Akiskan akiskan;
        akiskan.Hiz = 10.0;
        akiskan.Yogunluk = 1.2;

        double basinc = HesaplaBasinc(akiskan);

        std::cout << "Akiskanin Basinci: " << basinc << std::endl;
    }
    else if (cho == 2) {
            Material steel(50.2, 0.45);
            Part part1(20.0, 2.0);
            Part part2(20.0, 3.0);
            double heatTransfer1 = calculateHeatTransfer(steel.thermalConductivity, steel.specificHeat, part1.initialTemperature, 100.0, part1.length);
            double heatTransfer2 = calculateHeatTransfer(steel.thermalConductivity, steel.specificHeat, part2.initialTemperature, 100.0, part2.length);
            double totalHeatTransfer = heatTransfer1 + heatTransfer2;
            std::cout << heatTransfer1;
            std::cout << heatTransfer2;
            std::cout << totalHeatTransfer;
    }
    else if (cho == 3) {
        Vehicle car;

        double totalTime = 0.0;
        double timeStep = 0.1; 

        for (int i = 0; i < 10; ++i) { 
            car.acceleration = 2.0; 
            car.update(timeStep);

            std::cout << "Zaman: " << totalTime << " saniye, Pozisyon: " << car.position << " metre, Hız: " << car.velocity << " m/s" << std::endl;
            totalTime += timeStep;
        }
    }
    else if (cho == 4) {
        const double g = 9.81; 
        double h = 100.0;      
        double t = 0.0;        
        double dt = 0.01;      

        while (h > 0) {
            double v = sqrt(2 * g * h); 
            std::cout << "Zaman: " << t << " s, Yükseklik: " << h << " m" << std::endl;
            t += dt;
            h -= 0.5 * g * pow(dt, 2);
        }
        std::cout << "Cisim yere düştü!" << std::endl;

    }
    else if (cho == 5) {
        DCMotorSimulation motor(0.1, 0.01, 0.02); 

        double simulationTime = 10.0; 
        double timeStep = 0.01;      

        for (double t = 0; t < simulationTime; t += timeStep) {
            motor.setTargetSpeed(10 * sin(t)); 
            motor.update(timeStep);

            std::cout << "Zaman: " << t << " s, Pozisyon: " << motor.getCurrentPosition() << " m" << std::endl;
        }
    }
    else if (cho == 6) {
        double E = 2.1e11;  
        double L = 2.0;     
        double F = 1000.0;  
        double deformation = (F * L) / E;
        std::cout << "Uygulanan Kuvvet: " << F << " N" << std::endl;
        std::cout << "Çubuğun Deformasyonu: " << deformation << " metre" << std::endl;
    }
    else if (cho == 7) {
        double metalThickness = 1.0;
        double corrosionRate = 0.01; 
        double simulationYears = 10.0; 

        std::cout << "Başlangıçta Metal Kalınlığı: " << metalThickness << " mm" << std::endl;

        for (double t = 0; t < simulationYears; ++t) {
            metalThickness -= corrosionRate;

            std::cout << "Yıl " << t + 1 << " Sonunda Metal Kalınlığı: " << metalThickness << " mm" << std::endl;
        }
    }
    else if (cho == 8) {
        double pistonArea = 0.01; 
        double loadForce = 1000.0;  
        double initialPosition = 0.0; 
        double timeStep = 0.1;       
        double simulationTime = 10.0; 
        HydraulicSystem hydraulicSystem(pistonArea, loadForce, initialPosition);
        hydraulicSystem.simulate(timeStep, simulationTime);
    }
    else if (cho == 9) {
        double surfaceTensionCoefficient = 0.072; 
        SurfaceTensionSimulator simulator(surfaceTensionCoefficient);
        double dropletRadius = 0.001; 
        double surfaceTension = simulator.calculateSurfaceTension(dropletRadius);
        std::cout << "Sıvı Damla Yarıçapı: " << dropletRadius << " m" << std::endl;
        std::cout << "Yüzey Gerilimi: " << surfaceTension << " N/m" << std::endl;
    }
    else if (cho == 10) {
        double initialPollution = 50.0; 
        double emissionRate = 2.0;      
        double simulationTime = 24.0;   
        Environment environment(initialPollution);
        std::cout << "Başlangıç Kirlilik Seviyesi: " << initialPollution << " birim" << std::endl;
        std::cout << "Emisyon Hızı: " << emissionRate << " birim/saat" << std::endl;
        std::cout << "Simülasyon Süresi: " << simulationTime << " saat" << std::endl;
        environment.simulate(simulationTime, emissionRate);
    }
    return 0;
}


