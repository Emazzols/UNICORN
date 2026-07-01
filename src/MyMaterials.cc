#include "MyMaterials.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Exception.hh"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <locale>
#include <algorithm>
#include <utility>

using namespace CLHEP;

/*
 * Defines the physical meaning of the first column
 * in the input spectrum file.
 */
enum class SpectrumAxis
{
    Wavelength,   // first column = wavelength
    Energy        // first column = energy
};

namespace
{
    /*
     * Read an optical spectrum from file and convert
     * the X axis into Geant4 energies.
     *
     * Accepted file formats:
     *
     *   650.0 0.12
     *   650.0,0.12
     *   650.0 , 0.12
     *
     * The function automatically ignores:
     *
     * - headers
     * - comments
     * - arbitrary text
     * - invalid lines
     *
     * The output energies are automatically sorted
     * in strictly increasing order as required by Geant4.
     *
     * PARAMETERS
     * ----------
     *
     * filename:
     *   Path to the input file.
     *
     * energies:
     *   Output vector containing Geant4 energies.
     *
     * values:
     *   Output vector containing spectrum values
     *   associated to each energy.
     *
     * axisType:
     *   Specifies whether the first column contains:
     *
     *   - wavelengths
     *   - energies
     *
     * xUnit:
     *   Unit applied to the first column.
     *
     *   Examples:
     *
     *   - nm
     *   - mm
     *   - eV
     *   - keV
     *
     * valueUnit:
     *   Unit applied to the second column.
     */
    void ReadOpticalSpectrum(
        const std::string& filename,
        std::vector<G4double>& energies,
        std::vector<G4double>& values,
        SpectrumAxis axisType,
        G4double xUnit,
        G4double valueUnit = 1.0)
    {
        // Open input file
        std::ifstream fin(filename);

        if (!fin.is_open())
        {
            G4Exception(
                "ReadOpticalSpectrum",
                "FileNotFound",
                FatalException,
                ("Cannot open file: " + filename).c_str());
        }

        /*
         * Force standard C locale.
         *
         * This guarantees that decimal numbers
         * use '.' as decimal separator.
         */
        fin.imbue(std::locale::classic());

        /*
         * Temporary container storing:
         *
         * pair.first  -> energy
         * pair.second -> spectrum value
         *
         * Using pairs allows sorting energies
         * while preserving correspondence with
         * spectrum values.
         */
        std::vector<std::pair<G4double, G4double>> spectrum;

        std::string line;

        // Read file line by line
        while (std::getline(fin, line))
        {
            // Ignore empty lines
            if (line.empty())
                continue;

            /*
             * Replace commas with spaces.
             *
             * This allows formats such as:
             *
             *   500.0,1.58
             */
            std::replace(line.begin(), line.end(), ',', ' ');

            // Create parser stream
            std::stringstream ss(line);

            ss.imbue(std::locale::classic());

            G4double x;
            G4double val;

            /*
             * Attempt to read two numeric columns.
             *
             * Invalid or non-numeric lines
             * are automatically ignored.
             */
            if (!(ss >> x >> val))
                continue;

            // Apply unit to first column
            x *= xUnit;

            G4double energy;

            /*
             * Convert first column into energy.
             *
             * If the file contains wavelengths:
             *
             *     E = hc / lambda
             *
             * If the file already contains energies:
             *
             *     E = E
             */
            if (axisType == SpectrumAxis::Wavelength)
            {
                energy = (h_Planck * c_light) / x;
            }
            else
            {
                energy = x;
            }

            // Store energy-value pair
            spectrum.push_back(
            {
                energy,
                val * valueUnit
            });
        }

        // Reject files with no valid data
        if (spectrum.empty())
        {
            G4Exception(
                "ReadOpticalSpectrum",
                "EmptyFile",
                FatalException,
                ("No valid data found in file: " + filename).c_str());
        }

        /*
         * Sort spectrum by increasing energy.
         *
         * This guarantees compatibility with Geant4
         * independently from the ordering used
         * in the original file.
         */
        std::sort(
            spectrum.begin(),
            spectrum.end(),
            [](const auto& a, const auto& b)
            {
                return a.first < b.first;
            });

        /*
         * Split sorted pairs into output vectors.
         */
        for (const auto& point : spectrum)
        {
            energies.push_back(point.first);
            values.push_back(point.second);
        }
    }
}

MyMaterials::MyMaterials()
{
  // default constructor
  // Nanomaterial (fast scintillator) and matrices (polymeric embedding)

  N  = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  H  = new G4Element("Hydrogen", "H", z = 1 , a = 1.01 * g / mole);
  C  = new G4Element("Carbon"  , "C", z = 6 , a = 12.01 * g / mole);
  O  = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  Lu = new G4Element ("Lutetium",  "Lu",  z = 71., a = 174.97 * g / mole);
  Si = new G4Element ("Silicon", "Si", z = 14 , a = 28.09 * g / mole);
  Bi = new G4Element("Bismuth",  "Bi", z = 83., a = 208.980 *g/mole);
  Ge = new G4Element("Germanium",  "Ge", z = 32., a = 72.64 *g/mole);
  B  = new G4Element(" Boron", "B", z=5 , a=10.811*g/mole);
  Na = new G4Element(" Sodium", "Na", z=11 , a=22.98977*g/mole);
  K  = new G4Element(" Potassium", "K", z=19 , a=39.0983*g/mole);
  Al = new G4Element(" Aluminum", "Al", z=13 , a=26.981538*g/mole);
  Cd = new G4Element ("Cadmiun",  "Cd",  z = 48, a = 112.414 * g / mole);
  Se = new G4Element ("Selenium",  "Se",  z = 34., a = 78.971 * g / mole);
  S  = new G4Element ("Sulfur",  "S",  z = 16., a = 32.065 * g / mole);
  F = new G4Element("Fluorine", "F", z=9., a=18.9984*g/mole);
  Te = new G4Element("Tellurium",  "Te", z = 52., a = 127.60 *g/mole);
  Ga = new G4Element("Gallium", "Ga", z=31., a=69.723*g/mole);
  In = new G4Element("Indium", "In", z=49., a=114.818*g/mole);

  fNistMan = G4NistManager::Instance();

}

MyMaterials::~MyMaterials()
{}

//--------------------------//
// VACUUM                      //
//--------------------------//
G4Material * MyMaterials::Vacuum()
{
  G4double density;
  G4Material* mat = new G4Material("Vacuum", density=1.0e-25*g/cm3, 1);
  mat->AddElement(H, 100.*perCent);

  AddMaterialToList(mat->GetName());

  return mat;
}


//--------------------------//
// AIR                      //
//--------------------------//
G4Material * MyMaterials::Air()
{
  G4double density;
  // Air - this is the external world
  G4Material* mat = new G4Material("Air", density=1.29*mg/cm3, 2);
  mat->AddElement(N, 70.*perCent);
  mat->AddElement(O, 30.*perCent);

  // NEW WAY OF ADDING PRE-DEFINED MATERIAL PROPERTIES
  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX", "Air");

  AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)

  return mat;
}
//--------------------------//


//--------------------------//
// AIR THIN LAYER           //
//--------------------------//
G4Material * MyMaterials::AirThinLayer()
{

  G4double density;
  // Air thin layer
  // we need another air element, same characteristics as the other one, but
  // different pointer, because we need to declare a surface between this air and esr
  G4Material* mat = new G4Material("AirThinLayer", density=1.29*mg/cm3, 2);
  mat->AddElement(N, 70.*perCent);
  mat->AddElement(O, 30.*perCent);

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", "Air");

  mat->SetMaterialPropertiesTable(myMPT2);
  AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)

  return mat;
}
//--------------------------//


//--------------------------//
// LYSO                     //
//--------------------------//
G4Material *MyMaterials::LYSO(G4double lightyield, G4double resolutionScale)
{
  G4double density;
  // LYSO

  // G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Material *mat = new G4Material ("LYSO", density = 7.10 * g / cm3, 3, kStateSolid);
  mat->AddElement(Lu, 2);
  mat->AddElement(Si, 1);
  mat->AddElement(O, 5);



  std::vector<G4double> RIND_Energy = {1.9076923077*eV, 1.9375*eV, 1.9682539683*eV, 2.0*eV, 2.0327868852*eV, 2.0666666667*eV, 2.1016949153*eV, 2.1379310345*eV, 2.1754385965*eV, 2.2142857143*eV, 2.2545454545*eV, 2.2962962963*eV, 2.3396226415*eV, 2.3846153846*eV, 2.431372549*eV, 2.48*eV, 2.5306122449*eV, 2.5833333333*eV, 2.6382978723*eV, 2.6956521739*eV, 2.7555555556*eV, 2.8181818182*eV, 2.8837209302*eV, 2.9523809524*eV, 3.0243902439*eV, 3.1*eV, 3.1794871795*eV, 3.2631578947*eV, 3.3513513514*eV, 3.4444444444*eV, 3.5428571429*eV};
  std::vector<G4double> RIND_INDEX = {1.8002667721, 1.8006743824, 1.8011074531, 1.8015683722, 1.8020598293, 1.8025848653, 1.8031469314, 1.8037499611, 1.804398457, 1.8050975982, 1.8058533729, 1.8066727422, 1.8075638479, 1.808536274, 1.8096013826, 1.810772748, 1.812066725, 1.8135032038, 1.815106624, 1.8169073619, 1.8189436598, 1.8212643626, 1.8239328873, 1.8270331239, 1.8306784679, 1.8350261081, 1.8403005296, 1.8468340056, 1.8551403831, 1.8660592056, 1.8810632555};
  std::vector<G4double> ABS_Energy = { 1.00 * eV , 2.82 * eV , 2.88 * eV , 2.95 * eV , 3.02 * eV  , 3.10 * eV  , 3.18 * eV  , 3.26 * eV , 4.08 * eV };
  std::vector<G4double> ABS_LENGTH = { 438.*mm   , 438.*mm   , 413.*mm   , 375.*mm   , 263.*mm    , 87.5 * mm  , 11.5 * mm  , 1.0 * mm  , 1.0 * mm  };

  std::vector<G4double> RayleighEner = {2.2543*eV, 4.08 *eV};
  std::vector<G4double> Rayleigh = {200*mm, 88*mm};



  G4MaterialPropertiesTable *MPT_LYSO = new G4MaterialPropertiesTable();


  std::vector<G4double> FAST_Energy    = {1.77169*eV, 1.77266*eV, 1.77558*eV, 1.77851*eV, 1.78145*eV, 1.78539*eV, 1.79033*eV, 1.7963*eV, 1.80231*eV, 1.80836*eV, 1.81445*eV, 1.82058*eV, 1.82882*eV, 1.83401*eV, 1.84553*eV, 1.85293*eV, 1.86147*eV, 1.869*eV, 1.87769*eV, 1.89308*eV, 1.90536*eV, 1.92007*eV, 1.93039*eV, 1.94901*eV, 1.95846*eV, 1.9668*eV, 1.97884*eV, 1.99102*eV, 2.00088*eV, 2.01209*eV, 2.02596*eV, 2.03617*eV, 2.04519*eV, 2.0569*eV, 2.06611*eV, 2.0794*eV, 2.09151*eV, 2.10239*eV, 2.112*eV, 2.1231*eV, 2.13431*eV, 2.14565*eV, 2.15566*eV, 2.16868*eV, 2.18038*eV, 2.19519*eV, 2.21171*eV, 2.2193*eV, 2.23464*eV, 2.23619*eV, 2.24395*eV, 2.25806*eV, 2.27234*eV, 2.28358*eV, 2.29493*eV, 2.30475*eV, 2.31631*eV, 2.32463*eV, 2.33134*eV, 2.33809*eV, 2.34487*eV, 2.35856*eV, 2.36719*eV, 2.37939*eV, 2.38642*eV, 2.40238*eV, 2.41134*eV, 2.424*eV, 2.43312*eV, 2.44047*eV, 2.44786*eV, 2.46278*eV, 2.47788*eV, 2.48741*eV, 2.49317*eV, 2.49702*eV, 2.50282*eV, 2.50865*eV, 2.5145*eV, 2.52038*eV, 2.52432*eV, 2.53223*eV, 2.5362*eV, 2.54619*eV, 2.55424*eV, 2.56031*eV, 2.56437*eV, 2.57049*eV, 2.57663*eV, 2.58487*eV, 2.59317*eV, 2.59734*eV, 2.60571*eV, 2.61414*eV, 2.61414*eV, 2.61837*eV, 2.62262*eV, 2.62475*eV, 2.62902*eV, 2.63331*eV, 2.63545*eV, 2.63976*eV, 2.64191*eV, 2.64841*eV, 2.65493*eV, 2.6593*eV, 2.66149*eV, 2.66588*eV, 2.67914*eV, 2.67914*eV, 2.68136*eV, 2.68136*eV, 2.68359*eV, 2.68805*eV, 2.68805*eV, 2.68805*eV, 2.69477*eV, 2.69477*eV, 2.69702*eV, 2.70153*eV, 2.70605*eV, 2.71286*eV, 2.71742*eV, 2.71971*eV, 2.722*eV, 2.722*eV, 2.72429*eV, 2.72889*eV, 2.72889*eV, 2.73351*eV, 2.73814*eV, 2.74279*eV, 2.74512*eV, 2.74979*eV, 2.75213*eV, 2.75447*eV, 2.75682*eV, 2.75917*eV, 2.76389*eV, 2.76389*eV, 2.76626*eV, 2.76626*eV, 2.77338*eV, 2.77576*eV, 2.78533*eV, 2.79255*eV, 2.79738*eV, 2.80223*eV, 2.80466*eV, 2.80709*eV, 2.80953*eV, 2.80953*eV, 2.81934*eV, 2.8218*eV, 2.82673*eV, 2.83168*eV, 2.84164*eV, 2.84916*eV, 2.85419*eV, 2.8643*eV, 2.86684*eV, 2.87449*eV, 2.87705*eV, 2.87961*eV, 2.88475*eV, 2.88733*eV, 2.8925*eV, 2.89509*eV, 2.90028*eV, 2.90549*eV, 2.90811*eV, 2.91073*eV, 2.91335*eV, 2.91335*eV, 2.91335*eV, 2.91861*eV, 2.92125*eV, 2.92125*eV, 2.92389*eV, 2.92654*eV, 2.92654*eV, 2.92919*eV, 2.92919*eV, 2.93185*eV, 2.93451*eV, 2.93717*eV, 2.93985*eV, 2.94252*eV, 2.9452*eV, 2.94789*eV, 2.94789*eV, 2.94789*eV, 2.95058*eV, 2.95868*eV, 2.96411*eV, 2.96955*eV, 2.96955*eV, 2.97228*eV, 2.97228*eV, 2.97228*eV, 2.97502*eV, 2.97502*eV, 2.97776*eV, 2.9805*eV, 2.9805*eV, 2.9805*eV, 2.98601*eV, 2.99154*eV, 2.99431*eV, 2.99431*eV, 2.99431*eV, 2.99708*eV, 2.99708*eV, 3.00544*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.00824*eV, 3.01385*eV, 3.0223*eV, 3.02797*eV, 3.02797*eV, 3.03081*eV, 3.03081*eV, 3.03081*eV, 3.03365*eV, 3.0365*eV, 3.03935*eV, 3.04221*eV, 3.04795*eV, 3.04795*eV, 3.05083*eV, 3.05371*eV, 3.05949*eV, 3.06239*eV, 3.06529*eV, 3.06529*eV, 3.0682*eV, 3.0682*eV, 3.07112*eV, 3.07696*eV, 3.08283*eV, 3.09464*eV, 3.09464*eV, 3.0976*eV, 3.10653*eV, 3.11252*eV, 3.11852*eV, 3.12757*eV, 3.13668*eV, 3.14583*eV, 3.15813*eV, 3.16741*eV, 3.17675*eV, 3.20828*eV, 3.23719*eV, 3.26664*eV, 3.28656*eV, 3.31351*eV, 3.34783*eV, 3.38287*eV};

  std::vector<G4double> FAST_COMPONENT = {0.011691, 0.011691, 0.011691, 0.0146138, 0.0146138, 0.0146138, 0.011691, 0.011691, 0.00876827, 0.00876827, 0.00584551, 0.00584551, 0.00584551, 0.00292276, 0.00876827, 0.0146138, 0.0146138, 0.0146138, 0.0204593, 0.023382, 0.0263048, 0.0204593, 0.0204593, 0.023382, 0.0292276, 0.0321503, 0.0350731, 0.0379958, 0.0379958, 0.0379958, 0.0350731, 0.0379958, 0.0409186, 0.0438413, 0.0526096, 0.0584551, 0.0643006, 0.0730689, 0.0730689, 0.0818372, 0.0906054, 0.0964509, 0.0993737, 0.105219, 0.111065, 0.122756, 0.125678, 0.146138, 0.160752, 0.146138, 0.157829, 0.163674, 0.184134, 0.192902, 0.20167, 0.219207, 0.230898, 0.242589, 0.25428, 0.265971, 0.274739, 0.292276, 0.306889, 0.315658, 0.321503, 0.350731, 0.368267, 0.385804, 0.397495, 0.415031, 0.432568, 0.458873, 0.482255, 0.496868, 0.514405, 0.529019, 0.549478, 0.564092, 0.581628, 0.593319, 0.602088, 0.616701, 0.637161, 0.660543, 0.681002, 0.71023, 0.736534, 0.756994, 0.777453, 0.806681, 0.844676, 0.868058, 0.891441, 0.9119, 0.938205, 0.955741, 0.984969, 1.0142, 1.03173, 1.05511, 1.07557, 1.11649, 1.13695, 1.15741, 1.17495, 1.19248, 1.21002, 1.22756, 1.27432, 1.2977, 1.31524, 1.32985, 1.36785, 1.40292, 1.39415, 1.4, 1.41754, 1.44092, 1.47015, 1.48476, 1.50814, 1.5286, 1.54906, 1.56952, 1.58998, 1.61921, 1.63967, 1.66597, 1.68935, 1.71566, 1.73904, 1.76242, 1.77996, 1.80042, 1.8238, 1.83549, 1.8618, 1.85303, 1.87933, 1.91733, 1.89979, 1.92902, 1.95825, 1.98163, 2.01378, 2.03424, 2.0547, 2.07808, 2.09562, 2.11023, 2.12484, 2.13361, 2.15407, 2.15699, 2.15992, 2.16576, 2.16868, 2.16868, 2.16284, 2.15699, 2.14823, 2.13946, 2.12484, 2.11023, 2.08977, 2.06639, 2.04593, 2.02839, 2.01086, 1.98455, 1.96409, 1.94948, 1.93194, 1.91733, 1.90271, 1.87641, 1.86472, 1.8501, 1.83841, 1.82088, 1.79749, 1.77119, 1.75073, 1.73027, 1.70689, 1.68058, 1.65428, 1.6309, 1.60167, 1.57244, 1.55491, 1.53152, 1.50522, 1.47891, 1.45261, 1.43215, 1.362, 1.40877, 1.38831, 1.33862, 1.31232, 1.27432, 1.28601, 1.25678, 1.21587, 1.19541, 1.17203, 1.14864, 1.12234, 1.10772, 1.06096, 1.08434, 1.0142, 0.987891, 0.967432, 0.938205, 0.9119, 0.879749, 0.853445, 0.82714, 0.786221, 0.739457, 0.765762, 0.681002, 0.660543, 0.716075, 0.637161, 0.60501, 0.581628, 0.552401, 0.531942, 0.505637, 0.485177, 0.458873, 0.435491, 0.412109, 0.356576, 0.379958, 0.309812, 0.336117, 0.280585, 0.25428, 0.175365, 0.157829, 0.207516, 0.13737, 0.119833, 0.0993737, 0.0759916, 0.0613779, 0.0526096, 0.0350731, 0.0263048, 0.011691, 0.00876827, 0.00876827, 0.011691, 0.011691, 0.011691, 0.00876827, 0.011691};

  /* If you want to load optical spectra from external files, uncomment the lines below.
     The same function can also be used for:
      - emission spectra
      - refractive index
      - absorption length
      - any other optical property
    Input files may contain either:
      - energies
      - wavelengths
    Please refer to the ReadOpticalSpectrum function defined at the beginning of this file 
    for a detailed description of the supported formats and options.*/
    
  // std::vector<G4double> FAST_Energy; //RIND_Energy, ABS_Energy
  // std::vector<G4double> FAST_COMPONENT; //RIND_INDEX, ABS_LENGHT

  // ReadOpticalSpectrum(
  //     "FilePath/LYSOEmission.txt",
  //     FAST_Energy, //RIND_Energy, ABS_Energy
  //     FAST_COMPONENT, //or RIND_INDEX, ABS_LENGHT
  //     SpectrumAxis::Energy, //or Wavelength
  //     eV //energy/wavelength units you have in the file
  //   ); 

  MPT_LYSO->AddProperty ("SCINTILLATIONCOMPONENT1", FAST_Energy , FAST_COMPONENT);
  MPT_LYSO->AddConstProperty ("SCINTILLATIONRISETIME1", 0.05 * ns);
  MPT_LYSO->AddConstProperty ("SCINTILLATIONTIMECONSTANT1", 36 * ns);

  MPT_LYSO->AddProperty ("RINDEX",        RIND_Energy , RIND_INDEX);
  MPT_LYSO->AddProperty ("ABSLENGTH",     ABS_Energy  , ABS_LENGTH);
  MPT_LYSO->AddProperty ("RAYLEIGH",      RayleighEner, Rayleigh);
  MPT_LYSO->AddConstProperty ("SCINTILLATIONYIELD", lightyield / MeV);
  MPT_LYSO->AddConstProperty ("RESOLUTIONSCALE", resolutionScale);

  mat->SetMaterialPropertiesTable (MPT_LYSO);

  AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)
  return mat;
}
//--------------------------//

//--------------------------//
// BC422                      //
//--------------------------//
G4Material *MyMaterials::BC422(G4double lightyield,G4double resolutionScale)
{
  G4double density;
  // atoms
  G4Material *mat = new G4Material ("BC422", density = 1.032 * g / cm3, 2, kStateSolid);
  mat->AddElement (C, 9);
  mat->AddElement (H, 10);

  //r_index
  std::vector<G4double> RIND_Energy = {1.9076923077*eV,3.5428571429*eV};
  std::vector<G4double> RIND_INDEX = {1.58,1.58 };

  // abs lenght
  std::vector<G4double> ABS_Energy = { 1.00 * eV , 4.08 * eV };
  std::vector<G4double> ABS_LENGTH = { 1000*m , 1000*m  };

  // emission (fast component)
  std::vector<G4double> FAST_Energy    = {2.722613071*eV,2.7598585789*eV,2.811839843*eV,2.8603260326*eV,2.9276366204*eV,2.9707472154*eV,3.0359941038*eV,3.0900362232*eV,3.1420643341*eV,3.1917748731*eV,3.2388620928*eV,3.2657856948*eV,3.2902574948*eV,3.3048247505*eV,3.3209984442*eV,3.3418135367*eV,3.358352206*eV,3.3735300556*eV,3.4043010808*eV,3.4324789254*eV,3.4547196096*eV,3.4804930851*eV,3.509951804*eV,3.5399134539*eV,3.5720996187*eV};
  std::vector<G4double> FAST_COMPONENT = {3.3277870216,7.4875207987,13.8103161398,20.9650582363,32.1131447587,39.6006655574,52.0798668885,62.0632279534,71.5474209651,81.1980033278,90.5158069884,95.840266223,99.5008319468,99.3344425957,97.0049916805,92.1797004992,84.8585690516,77.7038269551,64.3926788686,53.2445923461,43.4276206323,34.608985025,22.6289517471,11.1480865225,3.4941763727};

  G4MaterialPropertiesTable *MTP_BC422 = new G4MaterialPropertiesTable();

  MTP_BC422->AddProperty ("SCINTILLATIONCOMPONENT1", FAST_Energy , FAST_COMPONENT);
  MTP_BC422->AddProperty ("RINDEX",        RIND_Energy , RIND_INDEX);
  MTP_BC422->AddProperty ("ABSLENGTH",     ABS_Energy  , ABS_LENGTH);
  MTP_BC422->AddConstProperty ("SCINTILLATIONYIELD", lightyield / MeV);
  MTP_BC422->AddConstProperty ("RESOLUTIONSCALE", resolutionScale);
  MTP_BC422->AddConstProperty ("SCINTILLATIONTIMECONSTANT1", 1.26 * ns);
  //MTP_BC422->AddConstProperty ("YIELDRATIO", 1.0);
  MTP_BC422->AddConstProperty ("SCINTILLATIONRISETIME1", 0.03 * ns);
  mat->SetMaterialPropertiesTable (MTP_BC422);
  AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)

  return mat;
}




//--------------------------//
// BGO                      //
//--------------------------//
G4Material *MyMaterials::BGO(G4double lightyield,G4double resolutionScale)
{
  G4double density;
  //BGO

  G4Material *mat = new G4Material ("BGO", density = 7.13 *g/cm3, 3);       //diversa rispetto a quello
  mat->AddElement(Bi,4);                                                       // del borosillicato
  mat->AddElement(Ge,3);
  mat->AddElement(O,12);

  // BGO
  std::vector<G4double> BGO_FAST_Energy =  {3.77040860233365*eV, 3.61585235867366*eV, 3.53140914873929*eV, 3.45082003941488*eV, 3.37382797640524*eV, 3.30019563221144*eV, 3.23959353408589*eV, 3.21011952704912*eV, 3.17164484535450*eV, 3.14338888200161*eV, 3.10648875072251*eV, 3.07937684280798*eV, 3.05273407926177*eV, 3.03522639841588*eV, 3.01791912298039*eV, 3.00080810535235*eV, 2.99232479188995*eV, 2.98388930803080*eV, 2.95886652908419*eV, 2.94241676857188*eV, 2.91808189363016*eV, 2.90208126436090*eV, 2.89414623386354*eV, 2.87840631146125*eV, 2.87060004547027*eV, 2.85511450968779*eV, 2.83979450139379*eV, 2.83219663242935*eV, 2.81712093376006*eV, 2.80220551360321*eV, 2.79480656012035*eV, 2.78012587810855*eV, 2.76559862083641*eV, 2.74408988238856*eV, 2.72291311954838*eV, 2.68833523989256*eV, 2.65462511529282*eV, 2.62174993085634*eV, 2.58335857867149*eV, 2.55221449223097*eV, 2.50985297957746*eV, 2.48627176540535*eV, 2.45741150672681*eV, 2.44041463626819*eV, 2.42921357075978*eV, 2.41811438564507*eV, 2.39621846230377*eV, 2.38541807149944*eV, 2.36410766328197*eV, 2.35883919770931*eV, 2.34317419881547*eV, 2.32771588912495*eV, 2.31246020475876*eV, 2.30240050181524*eV, 2.28747377956199*eV, 2.27273935387774*eV, 2.25819353259492*eV, 2.24383271746619*eV, 2.22965340119699*eV, 2.21102394892551*eV, 2.19270322692913*eV, 2.17022518952495*eV, 2.14385233614028*eV, 2.11811275965203*eV, 2.09298392087343*eV, 2.06441012881293*eV, 2.04053218738336*eV, 2.01336364231269*eV, 1.99064545944511*eV, 1.96114003938097*eV, 1.93249680644733*eV, 1.89445158995505*eV, 1.86442015991700*eV, 1.83214918680001*eV, 1.79791749561215*eV, 1.75324804790486*eV, 1.71629203956926*eV, 1.67554054054054*eV, 1.65458786830393*eV, 1.64225165562914*eV};
  std::vector<G4double> BGO_FAST_COMPONENT = {0.002041054182111, 0.002219269223948, 0.009306690784171, 0.017049024956026, 0.035679677835476, 0.058594494072484, 0.085689649997522, 0.115419170773643, 0.144187811114547, 0.182030090942073, 0.216438371499533, 0.253915993265593, 0.292607006633558, 0.318039236557389, 0.350360039139377, 0.376922558513944, 0.388774626135172, 0.41704676738375, 0.457262675199632, 0.495893419979575, 0.551893642684025, 0.592447977140537, 0.630331781138524, 0.661748455280689, 0.686479417669235, 0.715303019136405, 0.744618692337455, 0.773940296285538, 0.804171678771895, 0.834926126581869, 0.857915322012682, 0.89327832355784, 0.929231443516381, 0.96697544598151, 1.00555566103964, 1.05165135618158, 1.09911711197258, 1.13107580839577, 1.16059577034025, 1.18463907821213, 1.21575319428376, 1.22014873336624, 1.22496661876872, 1.22260134565886, 1.21423414140848, 1.20061004996834, 1.18222520202002, 1.16237489297545, 1.15227568459528, 1.12614005973959, 1.10427044496425, 1.08151597124141, 1.05242461864591, 1.0178691516744, 0.98683474211791, 0.954749645526535, 0.921569354792846, 0.887316438523992, 0.857782078713543, 0.818897621360131, 0.754203111068263, 0.702159272294626, 0.650115831695716, 0.607805349841931, 0.549647113475471, 0.483288487364082, 0.431968931438889, 0.386442288077889, 0.344079370268331, 0.294156892073653, 0.233039005628303, 0.177815311598429, 0.158572967031361, 0.120988788026388, 0.080758881004827, 0.056636432724791, 0.049251447274109, 0.03100578584872, 0.031796033968499, 0};

  std::vector<G4double> BGO_RIND_Energy = {4.06524590163934*eV, 3.97403846153846*eV, 3.88805268109125*eV, 3.80454127032832*eV, 3.72566105769231*eV, 3.6489111241907*eV, 3.57629074127488*eV, 3.50551314673452*eV, 3.43843594009983*eV, 3.37295973884657*eV, 3.31081441922563*eV, 3.25006553079948*eV, 3.19232749742533*eV, 3.13581183611533*eV, 3.08202833706189*eV, 3.02931834839971*eV, 2.97909658817876*eV, 2.92982041587902*eV, 2.88281794931411*eV, 2.83665065202471*eV, 2.79256756756757*eV, 2.74922394678492*eV, 2.70779646210963*eV, 2.66702516670252*eV, 2.62802034760492*eV, 2.58959899749373*eV, 2.55281037677579*eV, 2.51654150598742*eV, 2.48178542834267*eV, 2.44749309119621*eV, 2.41460564751704*eV, 2.3821325648415*eV, 2.35096700796359*eV, 2.32017215568862*eV, 2.29059671162017*eV, 2.26135327375524*eV, 2.2332492795389*eV, 2.20544290288154*eV, 2.17870321560358*eV, 2.1522305155355*eV, 2.12675814751286*eV, 2.10152542372881*eV, 2.07723236723069*eV, 2.05315449577745*eV, 2.02996070726916*eV, 2.00696018128844*eV, 1.98479270049624*eV, 1.96279879689726*eV, 1.94159098026934*eV, 1.92053903345725*eV, 1.90022988505747*eV, 1.8800606520091*eV, 1.86059423769508*eV, 1.84125334125334*eV, 1.82257827429075*eV, 1.80401571366216*eV, 1.78608470181504*eV, 1.76825442099258*eV, 1.7510238666855*eV, 1.73388337295483*eV, 1.71731301939058*eV, 1.70082304526749*eV, 1.68487566245414*eV, 1.66899986539238*eV, 1.65364097092558*eV, 1.63834566596194*eV, 1.62354327615556*eV, 1.60879719735306*eV, 1.59452160493827*eV, 1.58029569207239*eV, 1.56651926721415*eV, 1.55278647463995*eV, 1.53948348646635*eV, 1.52621861152142*eV, 1.51336506774075*eV, 1.50054459639356*eV, 1.48811809889582*eV, 1.4757200666508*eV, 1.46369968126549*eV, 1.45170354759396*eV, 1.44006968641115*eV, 1.42845622119816*eV, 1.41719053606126*eV, 1.40594171674793*eV, 1.39502700270027*eV, 1.3841259209645*eV, 1.37354602858092*eV, 1.36297680554029*eV, 1.35271656120445*eV, 1.34246427024686*eV, 1.33250940354648*eV, 1.32256*eV, 1.31289707750953*eV, 1.30323733445449*eV, 1.2938536992591*eV, 1.28447114886564*eV, 1.27535486525406*eV, 1.26623774509804*eV, 1.25737754791603*eV, 1.24851475178733*eV, 1.2399*eV};
  std::vector<G4double> BGO_RIND_INDEX = {2.4104536947317, 2.3869564639829, 2.3660543649085, 2.3468322730031, 2.3295878519573, 2.3136068397993, 2.2991685736303, 2.2857011892325, 2.2734609527426, 2.2619804934999, 2.2514923950471, 2.2416081668034, 2.2325378914889, 2.223954043441, 2.2160460860945, 2.2085345585638, 2.2015903717376, 2.1949726074859, 2.1888356391057, 2.1829699111326, 2.1775151209435, 2.1722875602994, 2.1674139239205, 2.1627320447373, 2.1583570965174, 2.154145043028, 2.1502008387, 2.1463958426504, 2.1428259288252, 2.1393756324524, 2.1361327418796, 2.1329931495694, 2.1300374273656, 2.1271713204311, 2.1244689440351, 2.1218446410822, 2.1193667357903, 2.1169571140252, 2.1146788929522, 2.1124606086703, 2.110360684894, 2.1083135475364, 2.1063733796242, 2.104479845337, 2.1026832871231, 2.1009280461158, 2.0992609786396, 2.0976306167583, 2.0960806454585, 2.0945633649203, 2.0931195690328, 2.0917049566522, 2.0903576806548, 2.0890365144451, 2.0877771935845, 2.0865412805598, 2.0853622941703, 2.0842043336617, 2.0830988812732, 2.0820123492245, 2.0809743454594, 2.0799533960623, 2.0789773810997, 2.0780167628302, 2.0770978258298, 2.0761928095027, 2.0753265228672, 2.0744728397614, 2.0736552025019, 2.0728489909617, 2.0720763797395, 2.0713141389397, 2.070583265607, 2.0698618153911, 2.0691696900544, 2.0684861359458, 2.0678300347384, 2.0671817373677, 2.0665591742514, 2.0659437225705, 2.0653524245962, 2.0647676123467, 2.0642054978937, 2.0636493028824, 2.0631144624691, 2.0625850282712, 2.0620757075938, 2.0615713273612, 2.0610859122674, 2.0606050143678, 2.0601420175145, 2.0596831527673, 2.0592412017489, 2.0588030320555, 2.0583808588197, 2.0579621470167, 2.057558578408, 2.0571581791936, 2.0567721284896, 2.0563889802941, 2.0560194396164};
  G4MaterialPropertiesTable *MPT_BGO = new G4MaterialPropertiesTable();

  std::vector<G4double> ABS_Energy = { 1.00 * eV , 2.82 * eV , 2.88 * eV , 2.95 * eV , 3.02 * eV  , 3.10 * eV  , 3.18 * eV  , 3.26 * eV , 4.08 * eV };
  std::vector<G4double> ABS_LENGTH = { 438.*mm   , 438.*mm   , 413.*mm   , 375.*mm   , 263.*mm    , 87.5 * mm  , 11.5 * mm  , 1.0 * mm  , 1.0 * mm  };
  MPT_BGO->AddProperty ("SCINTILLATIONCOMPONENT1", BGO_FAST_Energy , BGO_FAST_COMPONENT);
  MPT_BGO->AddProperty ("RINDEX",        BGO_RIND_Energy , BGO_RIND_INDEX);
  MPT_BGO->AddProperty ("ABSLENGTH",     ABS_Energy  , ABS_LENGTH);
  //MPT_BGO->AddProperty ("RAYLEIGH",      RayleighEner, Rayleigh);
  MPT_BGO->AddConstProperty ("SCINTILLATIONYIELD", lightyield / MeV);
  MPT_BGO->AddConstProperty ("RESOLUTIONSCALE", resolutionScale);
  MPT_BGO->AddConstProperty ("SCINTILLATIONTIMECONSTANT1", 350 * ns);
  //MPT_BGO->AddConstProperty ("YIELDRATIO", 1.0);
  MPT_BGO->AddConstProperty ("SCINTILLATIONRISETIME1", 0.1 * ns);
  mat->SetMaterialPropertiesTable (MPT_BGO);
  AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)

  return mat;
}
//--------------------------//


//--------------------------//
// OPTICAL GREASE           //
//--------------------------//
G4Material *MyMaterials::OpticalGrease(G4double refr_index)
{
  G4double density;
  // Optical Grease
  G4Material *mat = new G4Material ("OpticalGrease", density = 1.0 * g / cm3, 3);
  mat->AddElement(C, 1);
  mat->AddElement(H, 1);
  mat->AddElement(O, 1);

  // Optical Grease
  std::vector<G4double> PhotonEnergy1 = { 0.0001 * eV, 2.820 * eV,  4.136 * eV };
  std::vector<G4double> RefractiveIndex1 = { refr_index,refr_index,refr_index };
  std::vector<G4double> Absorption1 = { 3.448*m, 3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m, 15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m, 45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m, 52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m, 30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m, 17.500*m, 14.500*m };//FIXME why 34 values???
  G4MaterialPropertiesTable *mt = new G4MaterialPropertiesTable();
  mt->AddProperty("RINDEX"   , PhotonEnergy1, RefractiveIndex1);
  //mt->AddProperty("ABSLENGTH", PhotonEnergy1, Absorption1     , nEntries1); WHAT IS THIS? 3 DOES NOT CORRESPOND TO 34
  mat->SetMaterialPropertiesTable (mt);
  AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)
  return mat;
}
//--------------------------//


//--------------------------//
// FUSED SILICA             //
//--------------------------//
G4Material *MyMaterials::FusedSilica()
{
  G4double density;
  G4Material* Fused_silica = new G4Material ("Fused_silica", density = 2.2 * g / cm3, 2);
  Fused_silica->AddElement (Si, 1);
  Fused_silica->AddElement (O, 2);

  // Fused silica

  std::vector<G4double> Energy = {1 *eV, 5 *eV };
  // std::vector<G4double> glass_ReR = {refr_index,refr_index,refr_index,refr_index,refr_index  };
  // std::vector<G4double> glass_ImR = {1./0.0000001,1./0.0000001,1./0.0000001,1./0.0000001,1./0.0000001};

  // std::vector<G4double> Energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };
  std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm};
  G4MaterialPropertiesTable* glass_mt = new G4MaterialPropertiesTable(); //why is here glass_mt?
  //glass_mt->AddProperty("RINDEX"   , ephotoni, glass_ReR);
  glass_mt->AddProperty("RINDEX", "Fused Silica");
  glass_mt->AddProperty("ABSLENGTH", Energy, glass_AbsLength);
  Fused_silica->SetMaterialPropertiesTable(glass_mt);
  AddMaterialToList(Fused_silica->GetName());  // always save the material in the list (function will check if it is not already there)
  return Fused_silica;
}
//--------------------------//

//--------------------------//
// Glass            //
//--------------------------//
G4Material *MyMaterials::Glass()
{
  G4double density;
  G4Material* Glass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
  Glass->AddElement(C, 91.533 * perCent);
  Glass->AddElement(H, 8.467 * perCent);

  std::vector<G4double> Energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };

  std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm, 420. * cm };
  auto glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH", Energy, glass_AbsLength);
  glass_mt->AddProperty("RINDEX", "Fused Silica");
  Glass->SetMaterialPropertiesTable(glass_mt);
}
//--------------------------//

//--------------------------//
// Aluminum         //
//--------------------------//
G4Material *MyMaterials::Aluminum()
{
  // Aluminum
  G4double density;
  G4Material* Aluminum = new G4Material("Aluminum", z = 13., a = 26.98 * g / mole, density = 2.7 * g / cm3);

  return Aluminum;   //AGGIUNTA MIA
}
//--------------------------//


//--------------------------//
// SILICIUM                 //
//--------------------------//
G4Material *MyMaterials::Silicium()
{
  G4double density;
  G4Material *Silicium = new G4Material ("Silicium", density = 2.3290 * g / cm3, 1);
  Silicium->AddElement (Si, 1);
  //silicon
  std::vector<G4double> siliconRindexEnergy =
  {
    1.2407002545,    1.2532325803,    1.2660206679,    1.2790724273,    1.2923960984,    1.3060002679,    1.3198938878,    1.3340862952,    1.3485872332,    1.3634068731,    1.3785558383,    1.3940452298,    1.4098866528,    1.4260922466,    1.4426747145,    1.4596473582,    1.4770241125,    1.4948195837,    1.5130490909,    1.5317287093,    1.5508753181,    1.5705066513,    1.5906413519,    1.6112990318,    1.6325003349,    1.654267006	,    1.6766219655,    1.6995893897,    1.7231947979,    1.7474651472,    1.772428935	,    1.7981163109,    1.8245591978,    1.8517914246,    1.8798488705,    1.9087696223,    1.9385941477,    1.9693654833,    2.0011294428,    2.0339348435,    2.0678337575,    2.1028817873,    2.1391383698,    2.1766671132,    2.2155361688,    2.2558186446,    2.2975930639,    2.3409438764,    2.3859620279,    2.4327455971,    2.481400509	,    2.5320413357,    2.5847921969,    2.6397877755,    2.6971744663,    2.7571116767,    2.8197733057,    2.8853494291,    2.954048225	,    3.0260981817,    3.1017506363,    3.1812827039,    3.2650006698,    3.3532439311,    3.4463895958,    3.54485787	,    3.6491183956,    3.7596977409,    3.8771882953,    4.0022588855,    4.135667515	,    4.2782767397,    4.4310723375,    4.5951861278,    4.7719240558,    4.962801018
  };
  std::vector<G4double> siliconRindexValue =
  {
    3.57,     3.574,    3.578,    3.582,    3.587,    3.592,    3.597,    3.602,    3.608,    3.614,    3.62 ,    3.626,    3.632,    3.638,    3.644,    3.65 ,    3.656,    3.662,    3.668,    3.674,    3.681,    3.688,    3.696,    3.705,    3.714,    3.723,    3.732,    3.741,    3.751,    3.762,    3.774,    3.787,    3.8  ,    3.815,    3.83 ,    3.844,    3.861,    3.879,    3.895,    3.916,    3.939,    3.962,    3.986,    4.015,    4.044,    4.077,    4.11 ,    4.15 ,    4.192,    4.239,    4.293,    4.348,    4.416,    4.491,    4.577,    4.676,    4.793,    4.925,    5.091,    5.305,    5.587,    5.976,    6.548,    6.863,    6.014,    5.483,    5.293,    5.179,    5.102,    5.074,    5.055,    4.426,    3.052,    2.129,    1.8  ,    1.694
  };
  std::vector<G4double> siliconAbsLenghtValue =
  {
    1.56E-02,    1.26E-02,    1.04E-02,    8.77E-03,    7.46E-03,    6.37E-03,    5.46E-03,    4.76E-03,    4.17E-03,    3.68E-03,    3.27E-03,    2.92E-03,    2.61E-03,    2.31E-03,    2.08E-03,    1.87E-03,    1.69E-03,    1.55E-03,    1.41E-03,    1.29E-03,    1.18E-03,    1.08E-03,    9.90E-04,    9.09E-04,    8.40E-04,    7.69E-04,    7.04E-04,    6.49E-04,    6.02E-04,    5.65E-04,    5.26E-04,    4.88E-04,    4.52E-04,    4.20E-04,    3.88E-04,    3.56E-04,    3.29E-04,    3.06E-04,    2.84E-04,    2.62E-04,    2.42E-04,    2.23E-04,    2.05E-04,    1.88E-04,    1.73E-04,    1.56E-04,    1.42E-04,    1.27E-04,    1.14E-04,    1.03E-04,    9.01E-05,    7.87E-05,    6.76E-05,    5.81E-05,    4.76E-05,    3.92E-05,    3.22E-05,    2.55E-05,    2.00E-05,    1.48E-05,    1.05E-05,    6.67E-06,    3.41E-06,    1.43E-06,    9.80E-07,    9.62E-07,    9.17E-07,    8.55E-07,    7.81E-07,    6.94E-07,    5.78E-07,    4.46E-07,    4.24E-07,    4.59E-07,    5.08E-07,    5.43E-07
  };
  // const G4int nBinsSilicon = sizeof(siliconRindexEnergy)/sizeof(G4double);
  assert(sizeof(siliconRindexEnergy) == sizeof(siliconRindexValue));
  G4MaterialPropertiesTable *Silicon_mt = new G4MaterialPropertiesTable();
  std::vector<G4double> siliconAbsLenghtEnergy = {1*eV,2*eV,3*eV};
  //std::vector<G4double> siliconAbsLenghtValue  = {0.1*mm,0.1*mm,0.1*mm};//TO DO CROSS CHECK: DO WE NEED IT? WAS COMMENTED
  Silicon_mt->AddProperty ("RINDEX", siliconRindexEnergy , siliconRindexValue);
  //Silicon_mt->AddProperty ("ABSLENGTH" , siliconAbsLenghtEnergy , siliconAbsLenghtValue);
  Silicium->SetMaterialPropertiesTable(Silicon_mt);
  // SilicioMPPC->SetMaterialPropertiesTable(Silicon_mt);
  AddMaterialToList(Silicium->GetName());  // always save the material in the list (function will check if it is not already there)
  return Silicium;
}
//--------------------------//


// //--------------------------//
// // SILICIUM MPPC            //
// //--------------------------//
// G4Material *MyMaterials::SilicioMPPC()
// {
//   G4double density;
//   G4int nelements;
//   G4Material *SilicioMPPC = new G4Material ("SilicioMPPC", density = 2.3290 * g / cm3, 1);
//   SilicioMPPC->AddElement (Si, 1);
//   //silicon
//   G4double siliconRindexEnergy[76] =
//   {
//     1.2407002545,    1.2532325803,    1.2660206679,    1.2790724273,    1.2923960984,    1.3060002679,    1.3198938878,    1.3340862952,    1.3485872332,    1.3634068731,    1.3785558383,    1.3940452298,    1.4098866528,    1.4260922466,    1.4426747145,    1.4596473582,    1.4770241125,    1.4948195837,    1.5130490909,    1.5317287093,    1.5508753181,    1.5705066513,    1.5906413519,    1.6112990318,    1.6325003349,    1.654267006	,    1.6766219655,    1.6995893897,    1.7231947979,    1.7474651472,    1.772428935	,    1.7981163109,    1.8245591978,    1.8517914246,    1.8798488705,    1.9087696223,    1.9385941477,    1.9693654833,    2.0011294428,    2.0339348435,    2.0678337575,    2.1028817873,    2.1391383698,    2.1766671132,    2.2155361688,    2.2558186446,    2.2975930639,    2.3409438764,    2.3859620279,    2.4327455971,    2.481400509	,    2.5320413357,    2.5847921969,    2.6397877755,    2.6971744663,    2.7571116767,    2.8197733057,    2.8853494291,    2.954048225	,    3.0260981817,    3.1017506363,    3.1812827039,    3.2650006698,    3.3532439311,    3.4463895958,    3.54485787	,    3.6491183956,    3.7596977409,    3.8771882953,    4.0022588855,    4.135667515	,    4.2782767397,    4.4310723375,    4.5951861278,    4.7719240558,    4.962801018
//   };
//   G4double siliconRindexValue[76] =
//   {
//     3.57,     3.574,    3.578,    3.582,    3.587,    3.592,    3.597,    3.602,    3.608,    3.614,    3.62 ,    3.626,    3.632,    3.638,    3.644,    3.65 ,    3.656,    3.662,    3.668,    3.674,    3.681,    3.688,    3.696,    3.705,    3.714,    3.723,    3.732,    3.741,    3.751,    3.762,    3.774,    3.787,    3.8  ,    3.815,    3.83 ,    3.844,    3.861,    3.879,    3.895,    3.916,    3.939,    3.962,    3.986,    4.015,    4.044,    4.077,    4.11 ,    4.15 ,    4.192,    4.239,    4.293,    4.348,    4.416,    4.491,    4.577,    4.676,    4.793,    4.925,    5.091,    5.305,    5.587,    5.976,    6.548,    6.863,    6.014,    5.483,    5.293,    5.179,    5.102,    5.074,    5.055,    4.426,    3.052,    2.129,    1.8  ,    1.694
//   };
//   G4double siliconAbsLenghtValue[76] =
//   {
//     1.56E-02,    1.26E-02,    1.04E-02,    8.77E-03,    7.46E-03,    6.37E-03,    5.46E-03,    4.76E-03,    4.17E-03,    3.68E-03,    3.27E-03,    2.92E-03,    2.61E-03,    2.31E-03,    2.08E-03,    1.87E-03,    1.69E-03,    1.55E-03,    1.41E-03,    1.29E-03,    1.18E-03,    1.08E-03,    9.90E-04,    9.09E-04,    8.40E-04,    7.69E-04,    7.04E-04,    6.49E-04,    6.02E-04,    5.65E-04,    5.26E-04,    4.88E-04,    4.52E-04,    4.20E-04,    3.88E-04,    3.56E-04,    3.29E-04,    3.06E-04,    2.84E-04,    2.62E-04,    2.42E-04,    2.23E-04,    2.05E-04,    1.88E-04,    1.73E-04,    1.56E-04,    1.42E-04,    1.27E-04,    1.14E-04,    1.03E-04,    9.01E-05,    7.87E-05,    6.76E-05,    5.81E-05,    4.76E-05,    3.92E-05,    3.22E-05,    2.55E-05,    2.00E-05,    1.48E-05,    1.05E-05,    6.67E-06,    3.41E-06,    1.43E-06,    9.80E-07,    9.62E-07,    9.17E-07,    8.55E-07,    7.81E-07,    6.94E-07,    5.78E-07,    4.46E-07,    4.24E-07,    4.59E-07,    5.08E-07,    5.43E-07
//   };
//   const G4int nBinsSilicon = sizeof(siliconRindexEnergy)/sizeof(G4double);
//   assert(sizeof(siliconRindexEnergy) == sizeof(siliconRindexValue));
//   G4MaterialPropertiesTable *Silicon_mt = new G4MaterialPropertiesTable();
//   G4double siliconAbsLenghtEnergy[] = {1*eV,2*eV,3*eV};
//   //G4double siliconAbsLenghtValue[]  = {0.1*mm,0.1*mm,0.1*mm};
//   Silicon_mt->AddProperty ("RINDEX", siliconRindexEnergy , siliconRindexValue, nBinsSilicon );
//   Silicon_mt->AddProperty ("ABSLENGTH" , siliconAbsLenghtEnergy , siliconAbsLenghtValue , 3);
//   // Silicio->SetMaterialPropertiesTable(Silicon_mt);
//   SilicioMPPC->SetMaterialPropertiesTable(Silicon_mt);
//   AddMaterialToList(SilicioMPPC->GetName());  // always save the material in the list (function will check if it is not already there)
//   return SilicioMPPC;
// }
// //--------------------------//


//--------------------------//
// BOROSILICATE GLASS       //
//--------------------------//
G4Material *MyMaterials::BorosilicateGlass()
{
  G4double density;
  //Borosilicate glass

  G4Material *Borosilicate_glass = new G4Material ("Borosilicate glass", density = 2.23 * g / cm3, 6); // at 25 °C
  Borosilicate_glass->AddElement (Na, 8);
  Borosilicate_glass->AddElement (K, 8);
  Borosilicate_glass->AddElement (Al, 4);
  Borosilicate_glass->AddElement (B, 26);
  Borosilicate_glass->AddElement (Si, 81);
  Borosilicate_glass->AddElement (O, 202);

  // Borosilicate glass
  //*************** GEANT4 V11.1.2 DEFINITION
  std::vector<G4double> wavelength_boro =
  {
    300*nm,321*nm,342*nm,363*nm,384*nm,405*nm,426*nm,447*nm,468*nm,489*nm,510*nm,531*nm,552*nm,573*nm,594*nm,615*nm,636*nm,657*nm,678*nm,699*nm,720*nm,741*nm,762*nm,783*nm,804*nm
  };
  std::vector<G4double> energy_boro(wavelength_boro.begin(), wavelength_boro.end());
  G4int length = wavelength_boro.size();
  for(int i=0; i<length; i++) { energy_boro.at(length-1-i) = 1239.8 / (wavelength_boro.at(i)/nm);}
                                    //we have to subtract 1, energy elements need to be in increasing
  for(int i=0; i<length; i++) { energy_boro.at(i) = energy_boro.at(i) * eV;}

  std::vector<G4double> RefrIndex_boro_reversed = {1.545242213,1.539122034,1.534071673,1.529892844,1.526410935,1.523484459,1.521002202,1.518877568,1.517043123,1.515446102,1.514044923,1.512806553,1.511704525,1.510717457,1.509827932,1.509021652,1.508286793,1.507613507,1.506993539,1.506419926,1.505886762,1.505389009,1.504922346,1.504483055,1.504067914};
  G4int length_refrIndex = RefrIndex_boro_reversed.size();
  std::vector<G4double> RefrIndex_boro(RefrIndex_boro_reversed.begin(), RefrIndex_boro_reversed.end());
  for(int i=0; i<length_refrIndex; i++) { RefrIndex_boro.at(length-1-i) = RefrIndex_boro_reversed.at(i);}
  // for(int i=0; i<length_refrIndex; i++) { G4cout << "E: " << energy_boro.at(i)/eV << " eV, n: " << RefrIndex_boro.at(i) << "." << G4endl;}

  G4MaterialPropertiesTable* Borosilicate = new G4MaterialPropertiesTable();

  Borosilicate->AddProperty("RINDEX", energy_boro, RefrIndex_boro);
  Borosilicate_glass->SetMaterialPropertiesTable(Borosilicate);
  AddMaterialToList(Borosilicate_glass->GetName());  // always save the material in the list (function will check if it is not already there)
  return Borosilicate_glass;
}
//--------------------------//


//--------------------------//
// EPOXY                    //
//--------------------------//
G4Material *MyMaterials::Epoxy(G4double refr_index)
{
  G4double density;
  // Epoxy
  G4Material *Epoxy = new G4Material ("Epoxy", density = 1.0 * g / cm3, 3);
  Epoxy->AddElement (C, 1);
  Epoxy->AddElement (H, 1);
  Epoxy->AddElement (O, 1);

  //Epoxy
  std::vector<G4double> RIND1_Energy    = { 1.0 * eV, 1.84 * eV, 4.08 * eV };
  std::vector<G4double> RIND2_INDEX     = { refr_index,refr_index,refr_index };
  G4MaterialPropertiesTable *epoxy_mt = new G4MaterialPropertiesTable();
  epoxy_mt->AddProperty("RINDEX", RIND1_Energy, RIND2_INDEX);
  Epoxy->SetMaterialPropertiesTable (epoxy_mt);
  AddMaterialToList(Epoxy->GetName());  // always save the material in the list (function will check if it is not already there)
  return Epoxy;
}
//--------------------------//


//--------------------------//
// TEFLON                   //
//--------------------------//
G4Material *MyMaterials::Teflon()
{
  G4double density;
  //Teflon
  G4Material* teflon = new G4Material("Teflon", density=2.2*g/cm3,2);
  teflon->AddElement(C,2);
  teflon->AddElement(F,4);
  AddMaterialToList(teflon->GetName());  // always save the material in the list (function will check if it is not already there)
  return teflon;
}
//--------------------------//

//--------------------------//
// CdSeCdS in PMMA          //
//--------------------------//
G4Material *MyMaterials::CdSeCdSinPMMA(G4double host_fractMass, //implemented by Vojtech
                                         G4double NP_fractMass,
                                         G4double lightyield,
                                         G4double risetime,
                                         G4double decaytime,
                                         G4double rayleigh_length,
                                         G4double mie_length)
{
  G4double PMMA_density = 1.190 * g / cm3;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  elements.push_back("C");
  natoms.push_back(5);
  elements.push_back("H");
  natoms.push_back(8);
  elements.push_back("O");
  natoms.push_back(2);

  PMMA_density = 1.190 * g / cm3;

  fPMMA = fNistMan->ConstructNewMaterial("PMMA", elements, natoms, PMMA_density);

  elements.clear();
  natoms.clear();

  G4double CdSe_density = 5.816 * g / cm3;
  G4double CdS_density = 4.826 * g / cm3;
  G4double CdSeS_density = (CdSe_density+CdS_density)/2.0;
  G4double composite_density = 1/(((host_fractMass/PMMA_density)+(NP_fractMass/CdSeS_density))/100.0);


  elements.push_back("Cd");
  natoms.push_back(1);
  elements.push_back("Se");
  natoms.push_back(1);
  elements.push_back("S");
  natoms.push_back(1);

  G4Material* CdSeCdS = fNistMan->ConstructNewMaterial("CdSeCdS", elements, natoms, CdSeS_density);

  elements.clear();
  natoms.clear();
  //CdSe
  // G4Material *CdSeS = new G4Material("CdSeS", CdSeS_density, 3, kStateSolid);//is the density correct?
  // CdSeS->AddElement(Cd, 1);
  // CdSeS->AddElement(Se, 1);
  // CdSeS->AddElement(S, 1);


  //Nanocomposite
  G4Material *CdSeCdSinPMMA = new G4Material("CdSeCdSinPMMA", composite_density, 2);
  CdSeCdSinPMMA->AddMaterial(fPMMA, host_fractMass*perCent);
  CdSeCdSinPMMA->AddMaterial(CdSeCdS, NP_fractMass*perCent);

  // Spectra measured by Vojtech (CdSeCdS in PMMA 0.1% conc, 3040 um thickness)

  //********************************************* PL spectrum *********************************************
  std::vector<G4double> CdSeCdSinPMMA_emission_wavelength = {300,300.5,301,301.5,302,302.5,303,303.5,304,304.5,305,305.5,306,306.5,307,307.5,308,308.5,309,309.5,310,310.5,311,311.5,312,312.5,313,313.5,314,314.5,315,315.5,316,316.5,317,317.5,318,318.5,319,319.5,320,320.5,321,321.5,322,322.5,323,323.5,324,324.5,325,325.5,326,326.5,327,327.5,328,328.5,329,329.5,330,330.5,331,331.5,332,332.5,333,333.5,334,334.5,335,335.5,336,336.5,337,337.5,338,338.5,339,339.5,340,340.5,341,341.5,342,342.5,343,343.5,344,344.5,345,345.5,346,346.5,347,347.5,348,348.5,349,349.5,350,350.5,351,351.5,352,352.5,353,353.5,354,354.5,355,355.5,356,356.5,357,357.5,358,358.5,359,359.5,360,360.5,361,361.5,362,362.5,363,363.5,364,364.5,365,365.5,366,366.5,367,367.5,368,368.5,369,369.5,370,370.5,371,371.5,372,372.5,373,373.5,374,374.5,375,375.5,376,376.5,377,377.5,378,378.5,379,379.5,380,380.5,381,381.5,382,382.5,383,383.5,384,384.5,385,385.5,386,386.5,387,387.5,388,388.5,389,389.5,390,390.5,391,391.5,392,392.5,393,393.5,394,394.5,395,395.5,396,396.5,397,397.5,398,398.5,399,399.5,400,400.5,401,401.5,402,402.5,403,403.5,404,404.5,405,405.5,406,406.5,407,407.5,408,408.5,409,409.5,410,410.5,411,411.5,412,412.5,413,413.5,414,414.5,415,415.5,416,416.5,417,417.5,418,418.5,419,419.5,420,420.5,421,421.5,422,422.5,423,423.5,424,424.5,425,425.5,426,426.5,427,427.5,428,428.5,429,429.5,430,430.5,431,431.5,432,432.5,433,433.5,434,434.5,435,435.5,436,436.5,437,437.5,438,438.5,439,439.5,440,440.5,441,441.5,442,442.5,443,443.5,444,444.5,445,445.5,446,446.5,447,447.5,448,448.5,449,449.5,450,450.5,451,451.5,452,452.5,453,453.5,454,454.5,455,455.5,456,456.5,457,457.5,458,458.5,459,459.5,460,460.5,461,461.5,462,462.5,463,463.5,464,464.5,465,465.5,466,466.5,467,467.5,468,468.5,469,469.5,470,470.5,471,471.5,472,472.5,473,473.5,474,474.5,475,475.5,476,476.5,477,477.5,478,478.5,479,479.5,480,480.5,481,481.5,482,482.5,483,483.5,484,484.5,485,485.5,486,486.5,487,487.5,488,488.5,489,489.5,490,490.5,491,491.5,492,492.5,493,493.5,494,494.5,495,495.5,496,496.5,497,497.5,498,498.5,499,499.5,500,500.5,501,501.5,502,502.5,503,503.5,504,504.5,505,505.5,506,506.5,507,507.5,508,508.5,509,509.5,510,510.5,511,511.5,512,512.5,513,513.5,514,514.5,515,515.5,516,516.5,517,517.5,518,518.5,519,519.5,520,520.5,521,521.5,522,522.5,523,523.5,524,524.5,525,525.5,526,526.5,527,527.5,528,528.5,529,529.5,530,530.5,531,531.5,532,532.5,533,533.5,534,534.5,535,535.5,536,536.5,537,537.5,538,538.5,539,539.5,540,540.5,541,541.5,542,542.5,543,543.5,544,544.5,545,545.5,546,546.5,547,547.5,548,548.5,549,549.5,550,550.5,551,551.5,552,552.5,553,553.5,554,554.5,555,555.5,556,556.5,557,557.5,558,558.5,559,559.5,560,560.5,561,561.5,562,562.5,563,563.5,564,564.5,565,565.5,566,566.5,567,567.5,568,568.5,569,569.5,570,570.5,571,571.5,572,572.5,573,573.5,574,574.5,575,575.5,576,576.5,577,577.5,578,578.5,579,579.5,580,580.5,581,581.5,582,582.5,583,583.5,584,584.5,585,585.5,586,586.5,587,587.5,588,588.5,589,589.5,590,590.5,591,591.5,592,592.5,593,593.5,594,594.5,595,595.5,596,596.5,597,597.5,598,598.5,599,599.5,600,600.5,601,601.5,602,602.5,603,603.5,604,604.5,605,605.5,606,606.5,607,607.5,608,608.5,609,609.5,610,610.5,611,611.5,612,612.5,613,613.5,614,614.5,615,615.5,616,616.5,617,617.5,618,618.5,619,619.5,620,620.5,621,621.5,622,622.5,623,623.5,624,624.5,625,625.5,626,626.5,627,627.5,628,628.5,629,629.5,630,630.5,631,631.5,632,632.5,633,633.5,634,634.5,635,635.5,636,636.5,637,637.5,638,638.5,639,639.5,640,640.5,641,641.5,642,642.5,643,643.5,644,644.5,645,645.5,646,646.5,647,647.5,648,648.5,649,649.5,650,650.5,651,651.5,652,652.5,653,653.5,654,654.5,655,655.5,656,656.5,657,657.5,658,658.5,659,659.5,660,660.5,661,661.5,662,662.5,663,663.5,664,664.5,665,665.5,666,666.5,667,667.5,668,668.5,669,669.5,670,670.5,671,671.5,672,672.5,673,673.5,674,674.5,675,675.5,676,676.5,677,677.5,678,678.5,679,679.5,680,680.5,681,681.5,682,682.5,683,683.5,684,684.5,685,685.5,686,686.5,687,687.5,688,688.5,689,689.5,690,690.5,691,691.5,692,692.5,693,693.5,694,694.5,695,695.5,696,696.5,697,697.5};

  std::vector<G4double> CdSeCdSinPMMA_emission_energy(CdSeCdSinPMMA_emission_wavelength.begin(), CdSeCdSinPMMA_emission_wavelength.end());
  G4int length_PL = CdSeCdSinPMMA_emission_wavelength.size();
  for(int i=0; i<length_PL; i++) { CdSeCdSinPMMA_emission_energy.at(length_PL-1-i) = 1239.8 / CdSeCdSinPMMA_emission_wavelength.at(i) * eV;} //Note the index of the energy vector! Swaps the order of the energy vector -> need to do the same with PL emission values so that they match again!



  //PL emission before swapping
  std::vector<G4double> CdSeCdSinPMMA_emission_reversed = {0.019,0.019113636,0.019606061,0.021386364,0.021075758,0.019568182,0.014393939,0.011575758,0.012590909,0.016204545,0.017613636,0.01469697,0.014060606,0.014840909,0.015780303,0.014810606,0.014931818,0.015439394,0.015628788,0.016409091,0.018348485,0.019848485,0.019121212,0.01705303,0.015234848,0.014257576,0.016128788,0.018106061,0.018431818,0.016310606,0.014795455,0.015363636,0.016765152,0.016719697,0.014318182,0.011340909,0.009666667,0.009462121,0.01019697,0.012560606,0.014856061,0.015530303,0.014431818,0.014227273,0.014121212,0.014,0.011765152,0.011439394,0.011598485,0.014083333,0.014515152,0.015462121,0.014916667,0.017469697,0.017992424,0.019484848,0.017189394,0.014659091,0.010856061,0.011863636,0.014401515,0.01725,0.017295455,0.015712121,0.013962121,0.011727273,0.010818182,0.010833333,0.012409091,0.012909091,0.013984848,0.014689394,0.016825758,0.016560606,0.014,0.010674242,0.009833333,0.010287879,0.010606061,0.011219697,0.013151515,0.016030303,0.016272727,0.014931818,0.013416667,0.013287879,0.012765152,0.010787879,0.011431818,0.014893939,0.0175,0.015098485,0.012712121,0.012969697,0.014939394,0.015204545,0.015015152,0.015681818,0.01505303,0.015015152,0.015727273,0.017787879,0.017825758,0.016651515,0.015075758,0.014424242,0.013916667,0.012219697,0.010015152,0.007393939,0.007007576,0.007583333,0.009568182,0.011083333,0.012643939,0.013174242,0.013598485,0.013431818,0.012666667,0.012431818,0.013106061,0.012992424,0.011606061,0.010833333,0.013833333,0.017871212,0.019113636,0.015128788,0.012522727,0.012136364,0.014734848,0.015863636,0.015780303,0.014295455,0.013772727,0.01244697,0.011431818,0.011439394,0.015477273,0.019060606,0.019,0.016106061,0.015787879,0.015962121,0.014325758,0.011590909,0.011272727,0.013825758,0.017772727,0.020560606,0.018840909,0.015795455,0.013871212,0.014356061,0.012727273,0.011015152,0.011598485,0.015227273,0.017666667,0.018106061,0.016757576,0.0145,0.012522727,0.012666667,0.013393939,0.015734848,0.017090909,0.017575758,0.015886364,0.016295455,0.015393939,0.012727273,0.009712121,0.010386364,0.013409091,0.015659091,0.015143939,0.012598485,0.008613636,0.007015152,0.011310606,0.016651515,0.01755303,0.01380303,0.009674242,0.008227273,0.009265152,0.011181818,0.012787879,0.014386364,0.014901515,0.014583333,0.014128788,0.013659091,0.014477273,0.01630303,0.019757576,0.021780303,0.020977273,0.018454545,0.019227273,0.020924242,0.020068182,0.017636364,0.018075758,0.019340909,0.017628788,0.01230303,0.01055303,0.012765152,0.01555303,0.017219697,0.018219697,0.018954545,0.019765152,0.021,0.021212121,0.020393939,0.016545455,0.01394697,0.013704545,0.017772727,0.019143939,0.018765152,0.015765152,0.018681818,0.021409091,0.021022727,0.014333333,0.011,0.010598485,0.011954545,0.010121212,0.01,0.012265152,0.017121212,0.020030303,0.020674242,0.017090909,0.013136364,0.010818182,0.012401515,0.014530303,0.01605303,0.015833333,0.016174242,0.015742424,0.015560606,0.014242424,0.013272727,0.013204545,0.014704545,0.01880303,0.022939394,0.024840909,0.022121212,0.017742424,0.014189394,0.015090909,0.019371212,0.023984848,0.025378788,0.022916667,0.019772727,0.017393939,0.017068182,0.017757576,0.019022727,0.019734848,0.021136364,0.0205,0.019143939,0.01655303,0.016484848,0.017590909,0.018060606,0.017212121,0.017659091,0.020606061,0.021234848,0.019780303,0.016075758,0.016734848,0.0195,0.023234848,0.022007576,0.018962121,0.014931818,0.013719697,0.014765152,0.017234848,0.018136364,0.01830303,0.019575758,0.021515152,0.020969697,0.018477273,0.017901515,0.021272727,0.024560606,0.026416667,0.026181818,0.025083333,0.022681818,0.021560606,0.022159091,0.022537879,0.022030303,0.021416667,0.022310606,0.022439394,0.02175,0.022356061,0.025159091,0.027,0.026106061,0.025348485,0.027204545,0.030636364,0.031810606,0.029734848,0.025795455,0.02319697,0.022787879,0.024075758,0.024356061,0.022992424,0.021984848,0.023522727,0.026242424,0.024848485,0.020666667,0.016045455,0.016371212,0.018909091,0.023628788,0.027068182,0.0285,0.026598485,0.024106061,0.024181818,0.026621212,0.028621212,0.027424242,0.027189394,0.031272727,0.037651515,0.040515152,0.037742424,0.044871212,0.052030303,0.057257576,0.059075758,0.066992424,0.072969697,0.077060606,0.080757576,0.094477273,0.113636364,0.131825758,0.147166667,0.164757576,0.189522727,0.211265152,0.227431818,0.235234848,0.246954545,0.270856061,0.301977273,0.332287879,0.361818182,0.389121212,0.417537879,0.453106061,0.497583333,0.543598485,0.59194697,0.649257576,0.713378788,0.77594697,0.833575758,0.909681818,1.007674242,1.115393939,1.220583333,1.341242424,1.500727273,1.687606061,1.872189394,2.05725,2.257643939,2.475310606,2.715212121,2.995984848,3.328083333,3.685363636,4.029954545,4.353537879,4.713901515,5.128424242,5.602212121,6.095848485,6.570204545,7.033386364,7.485174242,7.962113636,8.444689394,8.920772727,9.348219697,9.777825758,10.21818182,10.66271212,11.01290152,11.26948485,11.47718182,11.67229545,11.81384848,11.86554545,11.79981061,11.62759091,11.3689697,11.04399242,10.70415909,10.33232576,9.962356061,9.568227273,9.160219697,8.725325758,8.301818182,7.867992424,7.427166667,6.980931818,6.565606061,6.202272727,5.883265152,5.59525,5.318734848,5.042462121,4.794431818,4.621901515,4.50780303,4.381878788,4.249689394,4.162189394,4.146727273,4.138742424,4.113424242,4.073234848,4.060128788,4.073340909,4.145121212,4.250931818,4.366825758,4.459325758,4.520590909,4.57430303,4.600780303,4.593659091,4.581318182,4.619689394,4.68994697,4.724409091,4.716674242,4.724462121,4.739409091,4.694106061,4.595060606,4.507871212,4.422280303,4.306227273,4.148371212,4.011409091,3.859295455,3.704363636,3.558757576,3.44830303,3.329075758,3.186227273,3.0355,2.880916667,2.727121212,2.571613636,2.451931818,2.331106061,2.218651515,2.103863636,1.987992424,1.864621212,1.740060606,1.649916667,1.580090909,1.526151515,1.470386364,1.418030303,1.358674242,1.296583333,1.23844697,1.177257576,1.11255303,1.047424242,1.002719697,0.977613636,0.955015152,0.924007576,0.879090909,0.837598485,0.808242424,0.790333333,0.771340909,0.746818182,0.708901515,0.665939394,0.636340909,0.630393939,0.633257576,0.623727273,0.598090909,0.57455303,0.555212121,0.535469697,0.507810606,0.484174242,0.467878788,0.464310606,0.459015152,0.452439394,0.443537879,0.436856061,0.429219697,0.419136364,0.402393939,0.385886364,0.373886364,0.361681818,0.337204545,0.312174242,0.299151515,0.308295455,0.313575758,0.303060606,0.279037879,0.259272727,0.245681818,0.243560606,0.245159091,0.243272727,0.238992424,0.236121212,0.238098485,0.230492424,0.209833333,0.186015152,0.171265152,0.16605303,0.169454545,0.171643939,0.167098485,0.157787879,0.148037879,0.139878788,0.129689394,0.12544697,0.124318182,0.127113636,0.12569697,0.125613636,0.119325758,0.106325758,0.090643939,0.081280303,0.080401515,0.078030303,0.074272727,0.06530303,0.060537879,0.058575758,0.059227273,0.057060606,0.052734848,0.048818182,0.045810606,0.043833333,0.038265152,0.034484848,0.031992424,0.034621212,0.036015152,0.038787879,0.039356061,0.04125,0.038492424,0.033530303,0.026613636,0.022651515,0.020590909,0.020037879,0.019939394,0.022257576,0.025856061,0.028037879,0.025939394,0.022931818,0.02025,0.020151515,0.019583333,0.019984848,0.019712121,0.019954545,0.019348485,0.018825758,0.018030303,0.017492424,0.016848485,0.017068182,0.018189394,0.018568182,0.016409091,0.013477273,0.013848485,0.016810606,0.019492424,0.019598485,0.018174242,0.015856061,0.014969697,0.015613636,0.016939394,0.016825758,0.015636364,0.015363636,0.01580303,0.016098485,0.015901515,0.015954545,0.016969697,0.01744697,0.016772727,0.015174242,0.015075758,0.015780303,0.015969697,0.013810606,0.011954545,0.011901515,0.013704545,0.015371212,0.015765152,0.015515152,0.015454545,0.016075758,0.01605303,0.015537879,0.014098485,0.012712121,0.011212121,0.011227273,0.013068182,0.015045455,0.014977273,0.013128788,0.012333333,0.013621212,0.014871212,0.014666667,0.013825758,0.014022727,0.013848485,0.013333333,0.012810606,0.013416667,0.013515152,0.013219697,0.013,0.013022727,0.012090909,0.011386364,0.012068182,0.013871212,0.014560606,0.014257576,0.013909091,0.014439394,0.014575758,0.014189394,0.013560606,0.013568182,0.013068182,0.012242424,0.011340909,0.011401515,0.011689394,0.012356061,0.013227273,0.014840909,0.016090909,0.01644697,0.015371212,0.013325758,0.011348485,0.011159091,0.012537879,0.013893939,0.014651515,0.015007576,0.015189394,0.014606061,0.013757576,0.013568182,0.012810606,0.011590909,0.011219697,0.012772727,0.014265152,0.014689394,0.014469697,0.015159091,0.016636364,0.017424242,0.01580303,0.012954545,0.010863636,0.011795455,0.014022727,0.016083333,0.015416667,0.012295455,0.010212121,0.011090909,0.014166667,0.015515152,0.015848485,0.017030303,0.018090909,0.016810606,0.0135,0.012060606,0.012598485,0.014613636,0.015219697,0.014757576,0.013083333,0.012659091,0.013121212,0.013787879,0.013174242,0.01319697,0.013522727,0.014121212,0.013916667,0.014795455,0.015909091,0.016272727,0.015674242,0.015189394,0.014462121,0.013431818,0.013590909,0.014621212,0.015151515,0.014833333,0.015575758,0.017083333,0.017409091,0.015659091,0.014015152,0.013022727,0.012060606,0.011287879,0.011931818,0.012840909,0.012742424,0.011462121,0.011568182,0.013537879,0.015969697,0.016416667,0.015356061,0.013939394,0.014,0.014280303,0.014325758,0.013515152,0.013136364,0.012969697,0.012727273,0.011780303,0.010537879,0.009515152,0.00930303,0.009689394,0.010393939,0.01130303,0.012492424,0.013075758,0.012734848,0.012015152};

  //PL emission after swapping in order to match it to the correct energy
  std::vector<G4double> CdSeCdSinPMMA_emission(CdSeCdSinPMMA_emission_reversed.begin(), CdSeCdSinPMMA_emission_reversed.end());
  for(int i=0; i<length_PL; i++) { CdSeCdSinPMMA_emission.at(length_PL-1-i) =  CdSeCdSinPMMA_emission_reversed.at(i);}
  // for(int i=0; i<length_PL; i++) { G4cout << CdSeCdSinPMMA_emission_energy.at(i)/eV << " eV,\t" << CdSeCdSinPMMA_emission.at(i) << "."<< G4endl; }



  //********************************************* Absorption spectrum *********************************************
  //Absorption coefficient acquired from transmission measurement: 23/08/24, Thomas bench, 0.1% CdSeCdS in PMMA, 3040 um thickness, TOTAL measuement (angular+diffuse)!
  std::vector<G4double> CdSeCdSinPMMA_absorption_wavelength =
  {
    800*nm,795*nm,790*nm,785*nm,780*nm,775*nm,770*nm,765*nm,760*nm,755*nm,750*nm,745*nm,740*nm,735*nm,730*nm,725*nm,720*nm,715*nm,710*nm,705*nm,700*nm,695*nm,690*nm,685*nm,680*nm,675*nm,670*nm,665*nm,660*nm,655*nm,650*nm,645*nm,640*nm,635*nm,630*nm,625*nm,620*nm,615*nm,610*nm,605*nm,600*nm,595*nm,590*nm,585*nm,580*nm,575*nm,570*nm,565*nm,560*nm,555*nm,550*nm,545*nm,540*nm,535*nm,530*nm,525*nm,520*nm,515*nm,510*nm,505*nm,500*nm,495*nm,490*nm,485*nm,480*nm,475*nm,470*nm,465*nm,460*nm,455*nm,450*nm,445*nm,440*nm,435*nm,430*nm,425*nm,420*nm,415*nm,410*nm,405*nm,400*nm
  };
  std::vector<G4double> CdSeCdSinPMMA_absorption_energy(CdSeCdSinPMMA_absorption_wavelength.begin(), CdSeCdSinPMMA_absorption_wavelength.end());
  G4int length_abs = CdSeCdSinPMMA_absorption_wavelength.size();
  for(int i=0; i<length_abs; i++) { CdSeCdSinPMMA_absorption_energy.at(i) = 1239.8 / (CdSeCdSinPMMA_absorption_wavelength.at(i)/nm);}
  for(int i=0; i<length_abs; i++) { CdSeCdSinPMMA_absorption_energy.at(i) = CdSeCdSinPMMA_absorption_energy.at(i) * eV;}
  std::vector<G4double> CdSeCdSinPMMA_absorption =
  {
    3.98253*cm,3.96437*cm,3.81906*cm,3.83596*cm,3.81346*cm,3.72597*cm,3.70999*cm,3.72597*cm,3.70999*cm,3.61681*cm,3.58178*cm,3.51841*cm,3.45253*cm,3.3624*cm,3.34057*cm,3.33191*cm,3.28085*cm,3.27666*cm,3.25175*cm,3.22314*cm,3.20299*cm,3.15172*cm,3.12475*cm,3.07582*cm,3.02832*cm,2.98921*cm,2.94762*cm,2.90379*cm,2.86118*cm,2.8292*cm,2.79479*cm,2.76722*cm,2.72535*cm,2.67324*cm,2.62023*cm,2.5985*cm,2.55344*cm,2.51488*cm,2.47494*cm,2.43617*cm,2.39853*cm,2.35522*cm,2.32205*cm,2.27706*cm,2.23569*cm,2.19375*cm,2.15134*cm,2.10135*cm,2.04827*cm,1.99927*cm,1.94766*cm,1.88502*cm,1.80524*cm,1.6808*cm,1.56372*cm,1.50761*cm,1.47375*cm,1.37404*cm,1.11852*cm,0.845933*cm,0.79804*cm,0.889162*cm,0.934892*cm,0.903578*cm,0.834417*cm,0.761701*cm,0.72794*cm,0.725376*cm,0.707822*cm,0.650009*cm,0.558336*cm,0.451463*cm,0.369896*cm,0.316159*cm,0.270152*cm,0.221741*cm,0.183908*cm,0.16984*cm,0.167452*cm,0.174393*cm,0.187622*cm
  };

  // for(int i=0; i<length_abs; i++) { G4cout << CdSeCdSinPMMA_absorption_energy.at(i)/eV << " eV,\t" << CdSeCdSinPMMA_absorption.at(i)/cm << " cm."<< G4endl; }

  //********************************************* Rayleigh scattering length *********************************************
  std::vector<G4double> CdSeCdSinPMMA_Rayleigh_wavelength = {800*nm,750*nm,700*nm, 650*nm, 600*nm, 550*nm, 500*nm, 450*nm, 400*nm};
  std::vector<G4double> CdSeCdSinPMMA_Rayleigh_energy(CdSeCdSinPMMA_Rayleigh_wavelength.begin(), CdSeCdSinPMMA_Rayleigh_wavelength.end());

  G4int length_rayleigh = CdSeCdSinPMMA_Rayleigh_wavelength.size();
  for(int i=0; i<length_rayleigh; i++) { CdSeCdSinPMMA_Rayleigh_energy.at(i) = 1239.8 / (CdSeCdSinPMMA_Rayleigh_wavelength.at(i)/nm);}
  for(int i=0; i<length_rayleigh; i++) { CdSeCdSinPMMA_Rayleigh_energy.at(i) = CdSeCdSinPMMA_Rayleigh_energy.at(i) * eV;}
  std::vector<G4double> CdSeCdSinPMMA_Rayleigh(CdSeCdSinPMMA_Rayleigh_wavelength.begin(), CdSeCdSinPMMA_Rayleigh_wavelength.end());
  G4double CdSeCdSinPMMA_Rayleigh_base = rayleigh_length*mm; // absorption length for scattering of 500 nm photons from config file
  // G4cout << "Rayleigh base for 500 nm: " << CdSeCdSinPMMA_Rayleigh_base/mm << " mm." << G4endl;
  for(int i=0; i<length_rayleigh; i++) { CdSeCdSinPMMA_Rayleigh.at(i) = CdSeCdSinPMMA_Rayleigh_base*pow((CdSeCdSinPMMA_Rayleigh_wavelength.at(i)/(500*nm)),4);}
  for(int i=0; i<length_rayleigh; i++) { CdSeCdSinPMMA_Rayleigh.at(i) = CdSeCdSinPMMA_Rayleigh.at(i) * mm;}
  // for(int i=0; i<length_rayleigh; i++) { G4cout << CdSeCdSinPMMA_Rayleigh_wavelength.at(i)/nm   << " nm,\t"
  //                                               << CdSeCdSinPMMA_Rayleigh_energy.at(i)/eV       << " eV,\t"
  //                                               << CdSeCdSinPMMA_Rayleigh.at(i)/mm              << "mm," <<G4endl; }

  //********************************************* Mie scattering length *********************************************
  std::vector<G4double> CdSeCdSinPMMA_Mie_energy = {1.5*eV,8*eV};
  std::vector<G4double> CdSeCdSinPMMA_Mie_length = {mie_length*mm,mie_length*mm};

  //********************************************* Filling Material properties table *********************************************
  G4MaterialPropertiesTable *MPT_CdSeCdSinPMMA = new G4MaterialPropertiesTable();
  MPT_CdSeCdSinPMMA->AddProperty("RAYLEIGH",  CdSeCdSinPMMA_Rayleigh_energy,  CdSeCdSinPMMA_Rayleigh);
  MPT_CdSeCdSinPMMA->AddProperty("RINDEX", "PMMA");
  MPT_CdSeCdSinPMMA->AddConstProperty("SCINTILLATIONYIELD", lightyield / MeV);
  MPT_CdSeCdSinPMMA->AddConstProperty("RESOLUTIONSCALE", 1);
  MPT_CdSeCdSinPMMA->AddConstProperty("SCINTILLATIONTIMECONSTANT1", decaytime * ns);
  MPT_CdSeCdSinPMMA->AddConstProperty("SCINTILLATIONRISETIME1", risetime * ns);
  MPT_CdSeCdSinPMMA->AddProperty("ABSLENGTH", CdSeCdSinPMMA_absorption_energy, CdSeCdSinPMMA_absorption);
  MPT_CdSeCdSinPMMA->AddProperty("SCINTILLATIONCOMPONENT1", CdSeCdSinPMMA_emission_energy, CdSeCdSinPMMA_emission);

  MPT_CdSeCdSinPMMA->AddProperty("MIEHG",CdSeCdSinPMMA_Mie_energy,CdSeCdSinPMMA_Mie_length);
  MPT_CdSeCdSinPMMA->AddConstProperty("MIEHG_FORWARD", 0.98);
  MPT_CdSeCdSinPMMA->AddConstProperty("MIEHG_BACKWARD", 0.9);
  MPT_CdSeCdSinPMMA->AddConstProperty("MIEHG_FORWARD_RATIO", 0.95);
  // MPT_CdSeCdSinPMMA->AddProperty("ABSLENGTH", energy, absWLSfiber);
  // MPT_CdSeCdSinPMMA->AddProperty("SCINTILLATIONCOMPONENT1", energy, emissionFib);
  // MPT_CdSeCdSinPMMA->AddProperty("WLSABSLENGTH", CdSeCdSinPMMA_absorption_energy,  CdSeCdSinPMMA_absorption);
  // MPT_CdSeCdSinPMMA->AddProperty("WLSCOMPONENT",  CdSeCdSinPMMA_emission_energy, CdSeCdSinPMMA_emission);
  // MPT_CdSeCdSinPMMA->AddConstProperty("WLSTIMECONSTANT",  decaytime * ns);
  CdSeCdSinPMMA->SetMaterialPropertiesTable(MPT_CdSeCdSinPMMA);  //AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)
  AddMaterialToList(CdSeCdSinPMMA->GetName());  // always save the material in the list (function will check if it is not already there)

  return CdSeCdSinPMMA;
}

//--------------------------//
// Cu-doped CdS  NCs        //
//--------------------------//
G4Material *MyMaterials::CuDopedCdSNCs(G4double lightyield,
                                       G4double risetime,
                                       G4double decaytime,
                                       G4double rayleigh_length,
                                       G4double mie_length)
{
  G4double density = 1.190 * g / cm3;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;


  elements.push_back("Cu");
  natoms.push_back(1);
  elements.push_back("Cd");
  natoms.push_back(1);
  elements.push_back("S");
  natoms.push_back(1);

  G4Material* CuDopedCdSNCs = fNistMan->ConstructNewMaterial("CuDopedCdSNCs", elements, natoms, density);

  elements.clear();
  natoms.clear();

  // Spectra acquired from email Sergio 11/09/24

  //********************************************* PL spectrum *********************************************
  std::vector<G4double> CuDopedCdSNCs_emission_wavelength =
  {
    400,401.0599976,401.9599915,403.0299988,403.9299927,405,406.0599976,406.9599915,408.0299988,408.9299927,410,411.0700073,412,413.0700073,414,415.0700073,416,417.0700073,418,419.0700073,420,421.0599976,421.9599915,423.0299988,423.9299927,425,426.0599976,426.9599915,428.0299988,428.9299927,430,431.0599976,431.9599915,433.0299988,433.9299927,435,436.0599976,436.9599915,438.0299988,438.9299927,440,441.0599976,441.9599915,443.0299988,443.9299927,445,446.0599976,446.9599915,448.0299988,448.9299927,450,451.0599976,451.9599915,453.0299988,453.9299927,455,456.0599976,456.9599915,458.0299988,458.9299927,460,461.0599976,461.9599915,463.0299988,463.9299927,465,466.0599976,466.9599915,468.0299988,468.9299927,470,471.0599976,471.9599915,473.0299988,473.9299927,475,476.0599976,476.9599915,478.0299988,478.9299927,480,481.0400085,481.9400024,482.980011,484.019989,485.0700073,485.9700012,487.0100098,488.0499878,488.9500122,490,491.0599976,491.9599915,493.0299988,493.9299927,495,496.0599976,496.9599915,498.0299988,498.9299927,500,501.0400085,501.9400024,502.980011,504.019989,505.0700073,505.9700012,507.0100098,508.0499878,508.9500122,510,511.0400085,511.9400024,512.9799805,514.0200195,515.0700073,515.9699707,517.0100098,518.0499878,518.9500122,520,521.039978,521.9400024,522.9799805,524.0200195,525.0700073,525.9699707,527.0100098,528.0499878,528.9500122,530,531.039978,531.9400024,532.9799805,534.0200195,535.0700073,535.9699707,537.0100098,538.0499878,538.9500122,540,541.0200195,542.0499878,542.9400024,543.9699707,545,546.0200195,547.0499878,547.9400024,548.9699707,550,551.039978,551.9400024,552.9799805,554.0200195,555.0700073,555.9699707,557.0100098,558.0499878,558.9500122,560,561.0200195,562.0499878,562.9400024,563.9699707,565,566.0200195,567.0499878,567.9400024,568.9699707,570,571.0200195,572.0499878,572.9400024,573.9699707,575,576.0200195,577.0499878,577.9400024,578.9699707,580,581.0200195,582.0499878,582.9400024,583.9699707,585,586.0200195,587.0499878,587.9400024,588.9699707,590,591.0200195,592.0499878,592.9400024,593.9699707,595,596.0200195,597.0499878,597.9400024,598.9699707,600,601.0200195,602.0499878,602.9400024,603.9699707,605,606.0200195,607.0499878,607.9400024,608.9699707,610,611.0100098,612.0200195,613.039978,614.0499878,615.0700073,615.9400024,616.9500122,617.9699707,618.9799805,620,621.0200195,622.0499878,622.9400024,623.9699707,625,626.0200195,627.0499878,627.9400024,628.9699707,630,631.0100098,632.0200195,633.039978,634.0499878,635.0700073,635.9400024,636.9500122,637.9699707,638.9799805,640,641,642,643,644,645,646,647,648,649,650,651.0100098,652.0200195,653.039978,654.0499878,655.0700073,655.9400024,656.9500122,657.9699707,658.9799805,660,661,662,663,664,665,666,667,668,669,670,671.0100098,672.0200195,673.039978,674.0499878,675.0700073,675.9400024,676.9500122,677.9699707,678.9799805,680,681,682,683,684,685,686,687,688,689,690,690.9799805,691.9699707,692.9500122,693.9400024,695.0700073,696.0499878,697.039978,698.0200195,699.0100098,700,701,702,703,704,705,706,707,708,709,710,710.9799805,711.9699707,712.9500122,713.9400024,715.0700073,716.0499878,717.039978,718.0200195,719.0100098,720,720.9799805,721.9699707,722.9500122,723.9400024,725.0700073,726.0499878,727.039978,728.0200195,729.0100098,730,730.9799805,731.9699707,732.9500122,733.9400024,735.0700073,736.0499878,737.039978,738.0200195,739.0100098,740,740.9699707,741.9400024,743.0499878,744.0200195,745,745.9699707,746.9400024,748.0499878,749.0200195,750,750.9799805,751.9699707,752.9500122,753.9400024,755.0700073,756.0499878,757.039978,758.0200195,759.0100098,760,760.9699707,761.9400024,763.0499878,764.0200195,765,765.9699707,766.9400024,768.0499878,769.0200195,770,770.9500122,772.0499878,773.0100098,773.9699707,775.0599976,776.0200195,776.9799805,777.9400024,779.039978,780,780.9699707,781.9400024,783.0499878,784.0200195,785,785.9699707,786.9400024,788.0499878,789.0200195,790,790.9500122,792.0499878,793.0100098,793.9699707,795.0599976,796.0200195,796.9799805,797.9400024,799.039978,800,800.9500122,802.0499878,803.0100098,803.9699707,805.0599976,806.0200195,806.9799805,807.9400024,809.039978,810,810.9400024,812.0200195,812.9699707,814.0499878,815,815.9400024,817.0200195,817.9699707,819.0499878,820,820.9400024,822.0200195,822.9699707,824.0499878,825,825.9400024,827.0200195,827.9699707,829.0499878,830,830.9400024,832.0200195,832.9699707,834.0499878,835,835.9400024,837.0200195,837.9699707,839.0499878,840,840.9400024,842.0200195,842.9699707,844.0499878,845,845.9400024,847.0200195,847.9699707,849.0499878,850
  };

  std::vector<G4double> CuDopedCdSNCs_emission_energy(CuDopedCdSNCs_emission_wavelength.begin(), CuDopedCdSNCs_emission_wavelength.end());
  G4int length_PL = CuDopedCdSNCs_emission_wavelength.size();
  for(int i=0; i<length_PL; i++) { CuDopedCdSNCs_emission_energy.at(length_PL-1-i) = 1239.8 / CuDopedCdSNCs_emission_wavelength.at(i) * eV;} //Note the index of the energy vector! Swaps the order of the energy vector -> need to do the same with PL emission values so that they match again!



  //PL emission before swapping
  std::vector<G4double> CuDopedCdSNCs_emission_reversed =
  {
    1.868326426,2.48844552,1.675605059,2.379570007,3.007236004,2.879714489,2.037539959,4.523386955,4.98180294,4.490431786,4.292017937,3.51782155,3.65536809,3.688799143,3.190287113,2.87420249,3.007236004,2.546918154,2.440072536,3.272601604,2.484923363,2.460353374,2.499296188,2.879714489,1.868326426,2.168498516,2.480193853,2.48844552,2.005225658,1.979718804,2.194015741,2.612567663,2.180490255,2.254148245,2.379570007,2.085392475,2.496575117,2.274251699,2.295774221,2.234824657,2.109426975,2.415238857,1.675605059,2.037539959,2.448853731,2.580034256,2.431828976,2.508009911,2.712652206,2.351047754,2.469016075,2.319831848,2.834985971,2.482619762,2.900141001,2.999890327,3.224040747,2.536501169,3.187365532,3.088555813,3.215903997,3.608418941,3.531903267,4.374304295,4.501073837,4.289769173,4.48458004,5.377121925,4.889665604,5.203336239,6.095004082,6.18032217,6.213265419,6.695167542,7.428163052,7.652744293,8.171180725,9.37639904,9.178685188,9.86615181,10.70029736,10.4881649,12.90405655,12.12986565,13.09384632,14.30244446,14.34796619,15.59538269,16.29011154,17.94819641,18.09921265,19.09013557,20.97824097,22.60079956,23.43188667,25.25296974,26.8716526,28.49222946,29.52752304,31.61343384,32.7875061,34.95040131,37.35085678,39.12216568,40.43054962,43.6749649,46.44608688,48.97306061,52.52813721,52.70277786,55.7042923,57.50831985,62.57281113,64.90190125,66.15106201,69.94702911,71.47058105,74.56881714,78.47618866,80.56371307,82.79179382,85.84270477,88.49617004,91.97481537,95.21169282,97.34546662,101.0234299,106.2957382,108.531662,112.4953308,119.1345291,120.744751,125.1162796,130.9517517,136.3162842,140.9816284,147.6646118,152.0889435,157.5614624,160.1523132,167.9828949,174.7745667,182.8384552,185.1908112,193.4621277,197.9829102,207.0571289,212.5282745,221.7150116,229.1262665,235.3797913,242.953598,248.5693054,258.7311707,264.6492615,272.6151428,278.9210815,289.236145,296.0498657,303.192688,312.8325195,320.1980286,329.4508362,338.6755676,354.1169434,358.6325378,364.9664917,370.0197144,382.397522,392.1913757,403.7982483,411.1602783,418.6937256,426.52005,432.697052,448.0336914,452.2511597,460.3809204,478.1431274,485.0139465,491.1722107,496.5923157,509.2787476,513.5109863,524.0141602,537.4335938,544.5465088,546.2584839,555.3515015,560.3180542,573.8560791,579.3674316,589.4493408,590.8861694,609.1478882,606.7459717,613.3156128,628.5527344,637.4498901,636.4197998,648.5053711,649.5574341,652.6795654,657.0599365,663.9658813,670.8930664,679.0429688,670.1090698,681.5250244,690.4559326,680.2566528,690.3988647,697.1241455,694.8913574,697.8970337,697.7017822,701.852478,702.7885132,698.4024048,696.4110718,698.4634399,704.5296631,708.8190308,697.5599365,699.2399902,696.0516357,706.3010864,693.2763062,693.3691406,689.0234985,694.1138916,683.071106,680.963562,682.5423584,678.114624,668.8065186,661.8127441,663.8015747,657.1782227,654.8486328,648.7861328,636.3806763,638.9143677,627.8238525,624.7485352,615.5124512,616.265564,604.7022095,598.5668945,589.3955688,587.1744385,584.7739258,573.4075928,560.967041,561.8007202,553.6105957,544.0111694,545.3566895,529.7822876,524.1898193,511.8441467,512.7609253,505.9934387,497.5263367,487.1916504,481.2566528,474.7557373,462.7705383,458.5465698,446.8354492,444.4602356,436.7011414,431.5666199,418.8241272,410.6746826,405.4238892,402.0358887,382.7636108,380.8435669,372.5072327,365.6131897,360.6204834,357.2059937,346.8800659,342.6121826,332.4187927,327.8023682,318.180603,314.8123779,301.1410828,300.6431274,289.4402161,283.6147156,278.2929077,268.7710571,269.0884399,261.4515991,262.7357483,253.6278076,251.3772583,245.7200165,242.5290985,237.9594116,230.5718842,227.9451904,223.7126312,220.9499359,221.0666351,214.7121124,211.8357544,202.2469788,200.5354309,200.6123657,196.9802704,189.2029724,190.2013245,185.7038574,180.2963409,180.423172,177.4539795,171.1910248,168.0210114,163.3243561,159.9327393,159.1349487,155.0574341,152.5959778,148.4855957,143.6813049,141.3370972,138.5527649,134.2216797,132.8881531,128.3503571,126.643959,122.5003357,121.2939911,116.7658691,114.6939316,109.8629837,108.5566406,107.8249512,104.4908142,101.345665,99.17029572,98.02076722,94.56712342,93.71589662,89.13721466,86.29449463,85.70561218,81.68987274,79.59046173,79.18454742,75.97911072,73.27114105,70.9524231,70.53179932,69.60944366,67.24902344,66.3196106,62.11545944,63.23698044,61.1125946,59.74396133,57.79744339,56.03620529,55.55830383,54.65289307,52.23947906,51.23920059,50.9994812,47.28290558,47.26118088,43.53568268,44.96585083,43.98997116,42.13246918,41.76716614,40.30423355,39.59818268,39.23284912,38.7863312,37.50775299,37.26822372,36.04843475,35.07011261,35.39810486,32.91128082,32.63080978,31.61722961,30.6207733,30.4963837,31.4497528,29.3874207,28.3410339,29.75175629,26.99609375,25.05507278,24.63080978,29.01242828,23.63418961,22.04019165,20.9986763,21.36751175,20.12571716,18.88822746,19.63176727,18.31938553,17.93077087,17.50576973,17.26994705,17.00659561,15.96106434,16.33841515,15.20619774,14.3113842,14.68634796,15.14833832,13.37578869,13.60460377,13.49791622,13.32212257,11.85343266,12.6014204,11.29477787,11.25824261,10.67938614,10.75000954,9.877477646,9.378833771,8.778608322,8.66305542,9.263046265,7.855166912,7.858525753,6.81993103,6.388833523,6.987400532,5.747304916,6.22273159,4.936470985,4.741261482,4.883392811,3.9985044,4.476061821,3.441959858,3.245271921,3.290190458,2.961504459,2.687860489
  };

  //PL emission after swapping in order to match it to the correct energy
  std::vector<G4double> CuDopedCdSNCs_emission(CuDopedCdSNCs_emission_reversed.begin(), CuDopedCdSNCs_emission_reversed.end());
  for(int i=0; i<length_PL; i++) { CuDopedCdSNCs_emission.at(length_PL-1-i) =  CuDopedCdSNCs_emission_reversed.at(i);}
  // for(int i=0; i<length_PL; i++) { G4cout << CuDopedCdSNCs_emission_energy.at(i)/eV << " eV,\t" << CuDopedCdSNCs_emission.at(i) << "."<< G4endl; }



  //********************************************* Absorption spectrum *********************************************
  //Absorption coefficient acquired from email Sergio 11/09/24
  std::vector<G4double> CuDopedCdSNCs_absorption_wavelength =
  {
    600.017395*nm,599.0145874*nm,598.0115356*nm,597.0081787*nm,596.0046387*nm,595.0007324*nm,593.996521*nm,592.9922485*nm,591.9875488*nm,590.982666*nm,590.0147705*nm,589.0093384*nm,588.0036621*nm,586.9978027*nm,585.9915161*nm,584.9851074*nm,584.015686*nm,583.0087891*nm,582.0015259*nm,580.9940796*nm,579.9863892*nm,579.0157471*nm,578.0075073*nm,576.9990845*nm,575.9902954*nm,574.9813843*nm,574.0095215*nm,573.000061*nm,571.9902954*nm,571.0177612*nm,570.0075684*nm,568.9971313*nm,567.9863281*nm,567.0128784*nm,566.0016479*nm,564.9902344*nm,564.0159912*nm,563.0040894*nm,561.9918213*nm,561.0169067*nm,560.0042725*nm,558.991333*nm,558.015686*nm,557.0022583*nm,555.9886475*nm,555.0122681*nm,553.9981689*nm,552.9838257*nm,552.0067749*nm,550.9920044*nm,550.0144653*nm,548.9992065*nm,547.9836426*nm,547.0054932*nm,545.9893799*nm,545.0108032*nm,543.9942627*nm,543.0151978*nm,541.99823*nm,541.0186768*nm,540.0013428*nm,538.9835205*nm,538.0033569*nm,536.9852295*nm,536.0045776*nm,534.9860229*nm,534.0048828*nm,532.9858398*nm,532.0042725*nm,530.9847412*nm,530.0028076*nm,528.9828491*nm,528.0003662*nm,527.0177612*nm,525.9971313*nm,525.0140991*nm,523.9929199*nm,523.0093994*nm,521.987854*nm,521.0039063*nm,519.9819336*nm,518.9974976*nm,518.0130005*nm,516.9902344*nm,516.005249*nm,514.9821167*nm,513.9967041*nm,513.0110474*nm,511.9872437*nm,511.0011597*nm,510.0148621*nm,508.990448*nm,508.0037231*nm,507.0168152*nm,505.9916687*nm,505.004364*nm,504.0168457*nm,502.9910583*nm,502.0031433*nm,501.0149536*nm,499.9885864*nm,499.0000305*nm,498.0111694*nm,496.9841614*nm,495.9949036*nm,495.0055542*nm,494.0159302*nm,492.9880676*nm,491.9980469*nm,491.0079041*nm,490.0174255*nm,488.988739*nm,487.9979553*nm,487.0068665*nm,486.015686*nm,484.986145*nm,483.9944763*nm,483.002655*nm,482.0106812*nm,481.0184326*nm,479.987915*nm,478.9952698*nm,478.0024414*nm,477.0094604*nm,476.0162964*nm,474.9847107*nm,473.9910889*nm,472.997345*nm,472.0033264*nm,471.0091858*nm,470.0148621*nm,468.9819946*nm,467.9872437*nm,466.9923096*nm,465.9971619*nm,465.0018616*nm,464.0063477*nm,463.0107117*nm,462.0147705*nm,461.0186768*nm,459.9841309*nm,458.9876709*nm,457.9910278*nm,456.9942322*nm,455.9971619*nm,454.9999695*nm,454.0026245*nm,453.0050049*nm,452.0072632*nm,451.0093384*nm,450.0112305*nm,449.0129089*nm,448.0144043*nm,447.0157471*nm,446.0168457*nm,445.0178528*nm,444.0186157*nm,443.0191956*nm,441.9812012*nm,440.9814758*nm,439.9814453*nm,438.9813538*nm,437.9810791*nm,437.0190735*nm,436.0183716*nm,435.0175476*nm,434.0165405*nm,433.0153809*nm,432.0140076*nm,431.0124817*nm,430.0107727*nm,429.0088806*nm,428.0067749*nm,427.0045776*nm,426.0021973*nm,424.9995728*nm,423.9967957*nm,422.9938965*nm,421.9907837*nm,420.9874878*nm,419.9840698*nm,419.019043*nm,418.0152283*nm,417.011261*nm,416.0071716*nm,415.0028687*nm,413.9984436*nm,412.9937744*nm,411.9889526*nm,410.9840088*nm,410.0175781*nm,409.0122681*nm,408.0067749*nm,407.0011902*nm,405.9953613*nm,404.9893799*nm,403.9832764*nm,403.0157471*nm,402.0092468*nm,401.002655*nm,399.9958801*nm,398.9889221*nm,397.981842*nm,397.0133362*nm,396.0058899*nm,394.9983521*nm,393.9906006*nm,392.9827271*nm,392.0134277*nm,391.005188*nm,389.9968262*nm,388.9882813*nm,388.0184021*nm,387.0095825*nm,386.0005798*nm,384.991394*nm,383.9820862*nm,383.0114441*nm,382.0018311*nm,380.9920654*nm,379.9821167*nm,379.0108948*nm,378.0006714*nm,376.9902954*nm,376.0186157*nm,375.0079041*nm,373.9970703*nm,372.9860535*nm,372.0137939*nm,371.0025635*nm,369.9910583*nm,369.0184021*nm,368.0066528*nm,366.9947815*nm,365.9827576*nm,365.009491*nm,363.9971313*nm,362.9846497*nm,362.0109863*nm,360.99823*nm,359.985321*nm,359.0112*nm,357.9980164*nm,356.9846497*nm,356.0101624*nm,354.996521*nm,353.9827576*nm,353.0078125*nm,351.9937439*nm,351.0185852*nm,350.0042419*nm,348.9897461*nm,348.0141296*nm,346.9993591*nm,345.9844666*nm,345.0084534*nm,343.9933167*nm,343.0170288*nm,342.0015564*nm,340.9859924*nm,340.0093384*nm,338.9934387*nm,338.0165405*nm,337.0003967*nm,335.9841309*nm,335.0068054*nm,333.9902649*nm,333.0126953*nm,331.9958801*nm,331.0180664*nm,330.0009766*nm,328.9837952*nm,328.0055542*nm,326.9880981*nm,326.0096436*nm,324.9919128*nm,324.0131836*nm,322.9952087*nm,322.0162354*nm,320.9980164*nm,320.0187683*nm,319.0002747*nm,317.9816284*nm,317.0020447*nm,315.9831848*nm,315.0033569*nm,313.9841919*nm,313.0041504*nm,311.9847717*nm,311.0044556*nm,309.9848022*nm,309.0042725*nm,307.984375*nm,307.0036011*nm,305.9834595*nm,305.0024414*nm,303.9820862*nm,303.000824*nm,302.0194702*nm,300.9987183*nm,300.0171204*nm,298.9961243*nm,298.0143127*nm,296.993103*nm,296.0110474*nm,294.9895935*nm,294.0072937*nm,292.9856262*nm,292.0031128*nm,290.9812012*nm,289.9984741*nm,289.015625*nm,287.9933472*nm,287.0102844*nm,285.987793*nm,285.0044861*nm,283.9817505*nm,282.9982605*nm,282.014679*nm
  };
  std::vector<G4double> CuDopedCdSNCs_absorption_energy(CuDopedCdSNCs_absorption_wavelength.begin(), CuDopedCdSNCs_absorption_wavelength.end());
  G4int length_abs = CuDopedCdSNCs_absorption_wavelength.size();
  for(int i=0; i<length_abs; i++) { CuDopedCdSNCs_absorption_energy.at(i) = 1239.8 / (CuDopedCdSNCs_absorption_wavelength.at(i)/nm); }
  for(int i=0; i<length_abs; i++) { CuDopedCdSNCs_absorption_energy.at(i) = CuDopedCdSNCs_absorption_energy.at(i) * eV; }
  std::vector<G4double> CuDopedCdSNCs_absorption =
  {
    165.825*cm,168.396*cm,166.588*cm,180.056*cm,171.319*cm,155.661*cm,163.576*cm,173.787*cm,166.78*cm,163.085*cm,166.652*cm,157.639*cm,159.257*cm,161.388*cm,162.353*cm,157.182*cm,164.505*cm,180.28*cm,183.169*cm,162.718*cm,166.908*cm,169.646*cm,159.316*cm,168.396*cm,175.685*cm,164.818*cm,155.941*cm,159.843*cm,172.957*cm,154.553*cm,156.109*cm,174.275*cm,153.298*cm,161.99*cm,169.184*cm,155.773*cm,162.779*cm,156.165*cm,166.396*cm,154.388*cm,160.138*cm,153.951*cm,160.969*cm,155.941*cm,165.131*cm,159.024*cm,151.164*cm,155.997*cm,150.171*cm,147.719*cm,148.833*cm,143.806*cm,142.252*cm,149.448*cm,144.476*cm,143.949*cm,139.869*cm,143.379*cm,140.685*cm,150.119*cm,146.178*cm,141.88*cm,141.142*cm,143.379*cm,141.834*cm,140.05*cm,137.959*cm,141.372*cm,145.103*cm,138.222*cm,139.286*cm,137.391*cm,144.092*cm,146.622*cm,144.668*cm,141.234*cm,144.813*cm,140.231*cm,142.298*cm,142.392*cm,138.974*cm,140.185*cm,143.711*cm,149.808*cm,138.708*cm,136.958*cm,132.447*cm,133.015*cm,134.166*cm,133.219*cm,131.286*cm,129.795*cm,134.665*cm,125.157*cm,131.128*cm,125.555*cm,124.12*cm,127.061*cm,117.919*cm,119.971*cm,115.259*cm,114.771*cm,110.423*cm,113.452*cm,110.564*cm,107.392*cm,101.518*cm,102.549*cm,99.5403*cm,104.448*cm,96.3386*cm,92.4227*cm,91.0471*cm,90.2711*cm,88.6496*cm,84.6083*cm,87.9139*cm,89.0495*cm,86.444*cm,84.2799*cm,81.573*cm,81.1159*cm,78.9626*cm,77.9143*cm,75.6347*cm,75.4114*cm,72.9661*cm,71.4065*cm,71.5241*cm,70.7205*cm,67.4789*cm,65.3074*cm,63.8199*cm,62.4345*cm,61.22*cm,58.0841*cm,57.024*cm,54.3684*cm,52.3625*cm,49.4528*cm,48.1*cm,46.5432*cm,44.167*cm,42.0299*cm,39.4706*cm,37.1891*cm,35.157*cm,33.1523*cm,30.8229*cm,28.4877*cm,25.9683*cm,23.6841*cm,21.5745*cm,19.6203*cm,17.7285*cm,15.8484*cm,14.2705*cm,12.8148*cm,11.36*cm,10.0802*cm,8.97525*cm,8.02377*cm,7.15018*cm,6.42562*cm,5.72155*cm,5.08506*cm,4.62301*cm,4.21195*cm,3.84223*cm,3.47291*cm,3.17778*cm,2.94365*cm,2.74261*cm,2.55853*cm,2.40439*cm,2.24166*cm,2.11749*cm,2.00651*cm,1.90171*cm,1.82131*cm,1.74424*cm,1.68072*cm,1.62856*cm,1.58264*cm,1.53764*cm,1.50064*cm,1.46918*cm,1.44296*cm,1.41989*cm,1.39881*cm,1.38223*cm,1.36691*cm,1.3533*cm,1.34167*cm,1.33164*cm,1.32325*cm,1.31511*cm,1.30885*cm,1.30306*cm,1.29487*cm,1.28977*cm,1.28379*cm,1.27754*cm,1.27171*cm,1.26416*cm,1.25669*cm,1.24896*cm,1.23964*cm,1.23088*cm,1.22187*cm,1.2108*cm,1.20047*cm,1.18765*cm,1.17625*cm,1.16328*cm,1.14904*cm,1.13614*cm,1.12218*cm,1.10668*cm,1.09244*cm,1.07776*cm,1.06263*cm,1.04741*cm,1.03402*cm,1.01921*cm,1.00437*cm,0.989716*cm,0.975025*cm,0.960387*cm,0.948419*cm,0.933575*cm,0.918873*cm,0.903489*cm,0.889785*cm,0.875338*cm,0.862565*cm,0.848971*cm,0.834673*cm,0.825241*cm,0.81492*cm,0.802641*cm,0.791006*cm,0.780307*cm,0.769421*cm,0.759884*cm,0.748803*cm,0.739368*cm,0.728787*cm,0.719056*cm,0.710724*cm,0.701382*cm,0.691967*cm,0.682192*cm,0.673214*cm,0.6652*cm,0.657016*cm,0.648582*cm,0.640217*cm,0.631312*cm,0.621312*cm,0.612201*cm,0.602643*cm,0.593736*cm,0.583568*cm,0.574389*cm,0.564539*cm,0.555687*cm,0.546897*cm,0.537112*cm,0.528156*cm,0.519162*cm,0.508981*cm,0.499034*cm,0.490282*cm,0.481578*cm,0.471411*cm,0.463509*cm,0.453476*cm,0.444098*cm,0.434898*cm,0.425976*cm,0.417221*cm,0.408185*cm,0.399275*cm,0.390622*cm,0.38293*cm,0.374844*cm,0.364278*cm,0.355505*cm,0.34797*cm,0.339715*cm,0.332474*cm,0.324466*cm,0.317474*cm,0.310475*cm,0.304188*cm,0.297765*cm,0.291427*cm,0.286325*cm,0.280749*cm,0.27541*cm,0.270829*cm,0.265947*cm,0.261769*cm,0.257463*cm,0.253365*cm,0.250474*cm,0.246179*cm,0.242541*cm,0.239534*cm,0.237242*cm,0.233915*cm,0.230726*cm,0.22875*cm,0.224837*cm,0.222501*cm,0.218711*cm,0.215984*cm,0.215534
  };
  G4double absScale = 1.0;
  for(int i=0; i<length_abs; i++) { CuDopedCdSNCs_absorption.at(i) = CuDopedCdSNCs_absorption.at(i)/absScale; }
  // for(int i=0; i<length_abs; i++) { G4cout << CuDopedCdSNCs_absorption_energy.at(i)/eV << " eV,\t" << CuDopedCdSNCs_absorption.at(i)/cm << " cm."<< G4endl; }

  //********************************************* Rayleigh scattering length *********************************************
  std::vector<G4double> CuDopedCdSNCs_Rayleigh_wavelength = {800*nm,750*nm,700*nm, 650*nm, 600*nm, 550*nm, 500*nm, 450*nm, 400*nm};
  std::vector<G4double> CuDopedCdSNCs_Rayleigh_energy(CuDopedCdSNCs_Rayleigh_wavelength.begin(), CuDopedCdSNCs_Rayleigh_wavelength.end());

  G4int length_rayleigh = CuDopedCdSNCs_Rayleigh_wavelength.size();
  for(int i=0; i<length_rayleigh; i++) { CuDopedCdSNCs_Rayleigh_energy.at(i) = 1239.8 / (CuDopedCdSNCs_Rayleigh_wavelength.at(i)/nm);}
  for(int i=0; i<length_rayleigh; i++) { CuDopedCdSNCs_Rayleigh_energy.at(i) = CuDopedCdSNCs_Rayleigh_energy.at(i) * eV;}
  std::vector<G4double> CuDopedCdSNCs_Rayleigh(CuDopedCdSNCs_Rayleigh_wavelength.begin(), CuDopedCdSNCs_Rayleigh_wavelength.end());
  G4double CuDopedCdSNCs_Rayleigh_base = rayleigh_length*mm; // absorption length for scattering of 500 nm photons from config file
  // G4cout << "Rayleigh base for 500 nm: " << CuDopedCdSNCs_Rayleigh_base/mm << " mm." << G4endl;
  for(int i=0; i<length_rayleigh; i++) { CuDopedCdSNCs_Rayleigh.at(i) = CuDopedCdSNCs_Rayleigh_base*pow((CuDopedCdSNCs_Rayleigh_wavelength.at(i)/(620*nm)),4);}
  for(int i=0; i<length_rayleigh; i++) { CuDopedCdSNCs_Rayleigh.at(i) = CuDopedCdSNCs_Rayleigh.at(i) * mm;}
  // for(int i=0; i<length_rayleigh; i++) { G4cout << CuDopedCdSNCs_Rayleigh_wavelength.at(i)/nm   << " nm,\t"
  //                                               << CuDopedCdSNCs_Rayleigh_energy.at(i)/eV       << " eV,\t"
  //                                               << CuDopedCdSNCs_Rayleigh.at(i)/mm              << "mm," <<G4endl; }

  //********************************************* Mie scattering length *********************************************
  std::vector<G4double> CuDopedCdSNCs_Mie_energy = {1.5*eV,8*eV};
  std::vector<G4double> CuDopedCdSNCs_Mie_length = {mie_length*mm,mie_length*mm};


  //********************************************* Filling Material properties table *********************************************
  G4MaterialPropertiesTable *MPT_CuDopedCdSNCs = new G4MaterialPropertiesTable();
  MPT_CuDopedCdSNCs->AddProperty("RAYLEIGH",  CuDopedCdSNCs_Rayleigh_energy,  CuDopedCdSNCs_Rayleigh);
  MPT_CuDopedCdSNCs->AddProperty("RINDEX", "PMMA");
  MPT_CuDopedCdSNCs->AddConstProperty("SCINTILLATIONYIELD", lightyield / MeV);
  MPT_CuDopedCdSNCs->AddConstProperty("RESOLUTIONSCALE", 1);
  MPT_CuDopedCdSNCs->AddConstProperty("SCINTILLATIONTIMECONSTANT1", decaytime * ns);
  MPT_CuDopedCdSNCs->AddConstProperty("SCINTILLATIONRISETIME1", risetime * ns);
  MPT_CuDopedCdSNCs->AddProperty("ABSLENGTH", CuDopedCdSNCs_absorption_energy, CuDopedCdSNCs_absorption);
  MPT_CuDopedCdSNCs->AddProperty("SCINTILLATIONCOMPONENT1", CuDopedCdSNCs_emission_energy, CuDopedCdSNCs_emission);

  MPT_CuDopedCdSNCs->AddProperty("MIEHG",CuDopedCdSNCs_Mie_energy,CuDopedCdSNCs_Mie_length);
  MPT_CuDopedCdSNCs->AddConstProperty("MIEHG_FORWARD", 0.98);
  MPT_CuDopedCdSNCs->AddConstProperty("MIEHG_BACKWARD", 0.9);
  MPT_CuDopedCdSNCs->AddConstProperty("MIEHG_FORWARD_RATIO", 0.95);
  CuDopedCdSNCs->SetMaterialPropertiesTable(MPT_CuDopedCdSNCs);  //AddMaterialToList(mat->GetName());  // always save the material in the list (function will check if it is not already there)
  AddMaterialToList(CuDopedCdSNCs->GetName());  // always save the material in the list (function will check if it is not already there)

  return CuDopedCdSNCs;
}


//--------------------------//

//--------------------------//
// PLASTIC                     //
//--------------------------//
G4Material *MyMaterials::Plastic()
{
  G4double density;
  G4Material* Plastic = new G4Material ("Plastic", density = 1.032 * g / cm3, 2);
  Plastic->AddElement (C, 9);
  Plastic->AddElement (H, 10);
  AddMaterialToList(Plastic->GetName());
  return Plastic;
}
//--------------------------//

//--------------------------//
// GALLIUM NITRIDE                     //
//--------------------------//
G4Material *MyMaterials::GalliumNitride()
{
  G4double density;
  G4Material* GalliumNitride = new G4Material ("GalliumNitride", density = 6.15 * g / cm3, 2);
  GalliumNitride->AddElement (Ga, 1);
  GalliumNitride->AddElement (N, 1);
  AddMaterialToList(GalliumNitride->GetName());
  return GalliumNitride;
}
//--------------------------//

//--------------------------//
// INDIUM GALLIUM NITRIDE                     //
//--------------------------//
G4Material *MyMaterials::InGaNitride()
{
  G4double density;
  G4Material* InGaNitride = new G4Material ("InGaNitride", density = 6.15  * g / cm3, 3);
  InGaNitride->AddElement (In, 3);  // 13% of gallium
  InGaNitride->AddElement (Ga, 17);
  InGaNitride->AddElement (N, 20);  // number of atoms are increased in order to be able to implement the 13% in integer numbers.
  AddMaterialToList(InGaNitride->GetName());
  return InGaNitride;
}
//--------------------------//

//--------------------------//
// ALUMINUM OXIDE                     //
//--------------------------//
G4Material *MyMaterials::AluminumOxide()
{
  G4double density;
  G4Material* AluminumOxide = new G4Material ("AluminumOxide", density = 3.95* g / cm3, 2);
  AluminumOxide->AddElement (Al, 2);
  AluminumOxide->AddElement (O, 3);
  AddMaterialToList(AluminumOxide->GetName());
  return AluminumOxide;
}


void MyMaterials::AddMaterialToList(std::string name)
{
  // add material to list, if it is not already there
  bool same = false;
  for(int j = 0 ; j < materials.size(); j++)
  {
    if(!materials[j].compare(name) ) same = true;
  }
  if(!same) materials.push_back(name);
}
