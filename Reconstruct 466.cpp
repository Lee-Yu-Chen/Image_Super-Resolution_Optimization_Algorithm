
// things to add : 
// - think of a method to do preprocessing of the normalisation - move the minus outside of the for loop
// - find out the reason and solution to the low CPU utilisation
// - change the logic flow -> for for for for if else, if possible
// - debugging for unsuccessful colours
// - getPixelRGB() can also be optimised : pixel.get() a large area and check pixel one by one, but not declare Quantum* pixel and pixel.get() every time
// - to make the program fast, do as much pre-processing as possible


// progress review meeting & pro forma
// meeting in study break
// **correlation as one milestone

// changes made:
// - develop a new method to access pixels
// - remove usage of function
// - data type long int -> long long int
// - float -> double
// 
// data overflowing problem




// **make it take 2 images with different sizes
// **make it more automatic

// 8 Dec 2023
// Demonstration of Cost Function.cpp
// Discussion of problem overcame while developing the program ( int, long int overflow problem)
// inspection of memory allocation for different data type as verification of the problem identified
// discussion of math theory behind the weighted sum/ML-like algorithm

// 11 Dec 2023
// 
// Demonstration of New Correlation.cpp
// dicussion of improvements to the program
//  - can be used for images with different sizes
//  - increase the level of automation in the program ( after collecting data in file, use a MATLAB script to plot it/use a C++ code to find the peaks)
// 
// the direction of the upcoming tasks : usage of C++ libraries to implement the mathematics of the ML-like algo
// the direction of the future implementation
// 
// 
// meeting in the study break
// progress review meeting
//
// 
// 
// 
// 15 December 2023
// 
// - Brief discussion of the project plan in study break
// - Update on current progress : in the middle of implementing the mathematical concept of the ML-like algo
// - Pro forma, project review meeting and meeting minutes
// 
// 
// 
// 
// 2-3 Jan 2024
// - Achieve image reconstruction with 466
// 
// 
// 
// 
// 
//




/*


#include <Magick++.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>

using namespace std;
using namespace Magick;
//using Eigen::MatrixXi;


void append(string txtfile, string text) {

    // reading the exising contents in the file
    fstream read(txtfile);

    // store the file contents in a variable
    string storage, buffer;
    while (getline(read, buffer)) {
        storage += buffer;
        storage += '\n';
    }

    // erase the last '\n'
    storage.erase(storage.length() - 1, 1);

    // add the new text at the end of the file
    storage += text;
    read.close();

    fstream write(txtfile);

    // write the existing contents + the added text to the file
    write << storage;
    write.close();
};


//void getPixelRGB(Image img, int x, int y, int* R, int* G, int* B) {
void getPixelRGB(Image img, int x, int y, Quantum* R, Quantum* G, Quantum* B) {
    Pixels pixel(img);

    Quantum* pix = pixel.get(x, y, 1, 1);

    *R = *pix;
    *G = *(pix + 1);
    *B = *(pix + 2);
};







int main(int argc, char** argv)
{
    InitializeMagick(*argv);

    int s = 2;      // scaling factor



    Image image;
    image.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");

    size_t w = image.size().width();
    size_t h = image.size().height();

    image.sample(Geometry(w / 5, h / 5));
    image.write("C:/Users/ACER/Desktop/asd.png");

    w = image.size().width();
    h = image.size().height();

    Image sampled = image;
    sampled.sample(Geometry(w / s, h / s));
    sampled.write("C:/Users/ACER/Desktop/sampled.png");



    /*
    // contructing image beyond the original resolution
    Image sampled;
    sampled.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");
    sampled.write("C:/Users/ACER/Desktop/asd.png");
    *





    // dimensions of the down-sampled image
    size_t width = sampled.size().width();
    size_t height = sampled.size().height();

    size_t width2 = width * s;
    size_t height2 = height * s;

    cout << "width : " << width2 << "\t height : " << height2 << endl << endl;


    // move the existing pixels on the down-sampled image to the reconstruction canva
    Image canva(Geometry(width2, height2), Color("#000000"));     // make an image with the original size
    Pixels samplePix(sampled);
    Pixels canvaPix(canva);
    Quantum* Qu1 = samplePix.get(0, 0, width, height);
    Quantum* Qu2 = canvaPix.get(0, 0, width2, height2);


    cout << "Copying the original pixels..." << endl;

    // 1st stage of building the reconstructed image, moving the existing pixels to the recon image 
    for (int i = 0; i <= height - 1; i++) {
        for (int j = 0; j <= width - 1; j++) {

            int I = ((i * width) + j) * 3;
            int I2 = ((i * 2 * width2) + j * 2) * 3;

            *(Qu2 + I2) = *(Qu1 + I);
            *(Qu2 + I2 + 1) = *(Qu1 + I + 1);
            *(Qu2 + I2 + 2) = *(Qu1 + I + 2);
        }
    }

    cout << "Done copying the original pixels" << endl;

    cout << "Syncing..." << endl;
    canvaPix.sync();

    cout << "Writing image..." << endl;
    canva.write("C:/Users/ACER/Desktop/canva.png");


    cout << "\nDeclaring 1st set of coefficients (4+4+4) ..." << endl;
    // the set of coefficients ai
    float R_a0 = 1;
    float R_a1 = -8.94488e-10;
    float R_a2 = 6.44377e-10;
    float R_a3 = -8.24457e-10;

    float G_a0 = 1;
    float G_a1 = 4.52474e-11;
    float G_a2 = -5.60703e-10;
    float G_a3 = 4.27008e-10;

    float B_a0 = 1;
    float B_a1 = -1.25965e-10;
    float B_a2 = -2.88537e-10;
    float B_a3 = -7.54881e-11;



    cout << "Declaring RGB images..." << endl;

    // testing recon
    Image redRecon(Geometry(width2, height2), Color("#000000"));
    Image greenRecon(Geometry(width2, height2), Color("#000000"));
    Image blueRecon(Geometry(width2, height2), Color("#000000"));
    Pixels redRecPix(redRecon);
    Pixels greenRecPix(greenRecon);
    Pixels blueRecPix(blueRecon);
    Quantum* Qur = redRecPix.get(0, 0, width2, height2);
    Quantum* Qug = greenRecPix.get(0, 0, width2, height2);
    Quantum* Qub = blueRecPix.get(0, 0, width2, height2);


    cout << "Generating 1st batch of reconstructed pixels" << endl;

    // 2nd stage of building the recon image, filling in the 1st batch reconstructed pixels
    for (int i = 1; i <= height2 - 3; i += 2) {
        for (int j = 1; j <= width2 - 3; j += 2) {          // i and j is the position of the interested pixel in canva.png

            int Ii = ((i * width2) + j) * 3;

            int I0 = (((i - 1) * width2) + j - 1) * 3;
            int I1 = (((i - 1) * width2) + j + 1) * 3;
            int I2 = (((i + 1) * width2) + j - 1) * 3;
            int I3 = (((i + 1) * width2) + j + 1) * 3;

            *(Qu2 + Ii) = (*(Qu2 + I0)) * R_a0 + (*(Qu2 + I1)) * R_a1 + (*(Qu2 + I2)) * R_a2 + (*(Qu2 + I3)) * R_a3;
            *(Qu2 + Ii + 1) = (*(Qu2 + I0 + 1)) * G_a0 + (*(Qu2 + I1 + 1)) * G_a1 + (*(Qu2 + I2 + 1)) * G_a2 + (*(Qu2 + I3 + 1)) * G_a3;
            *(Qu2 + Ii + 2) = (*(Qu2 + I0 + 2)) * B_a0 + (*(Qu2 + I1 + 2)) * B_a1 + (*(Qu2 + I2 + 2)) * B_a2 + (*(Qu2 + I3 + 2)) * B_a3;

            *(Qur + Ii) = (*(Qu2 + I0)) * R_a0 + (*(Qu2 + I1)) * R_a1 + (*(Qu2 + I2)) * R_a2 + (*(Qu2 + I3)) * R_a3;
            *(Qug + Ii + 1) = (*(Qu2 + I0 + 1)) * G_a0 + (*(Qu2 + I1 + 1)) * G_a1 + (*(Qu2 + I2 + 1)) * G_a2 + (*(Qu2 + I3 + 1)) * G_a3;
            *(Qub + Ii + 2) = (*(Qu2 + I0 + 2)) * B_a0 + (*(Qu2 + I1 + 2)) * B_a1 + (*(Qu2 + I2 + 2)) * B_a2 + (*(Qu2 + I3 + 2)) * B_a3;

        }
        cout << "1st :  " << i << endl;
    }

    cout << "Done generating the 1st batch of reconstructed pixels" << endl;

    cout << "Syncing..." << endl;
    canvaPix.sync();
    cout << "Writing image..." << endl;
    canva.write("C:/Users/ACER/Desktop/canva1.png");

    cout << "Syncing RGB..." << endl;
    //redRecPix.sync();
    //greenRecPix.sync();
    //blueRecPix.sync();
    //redRecon.write("C:/Users/ACER/Desktop/redRecon.png");
    //greenRecon.write("C:/Users/ACER/Desktop/greenRecon.png");
    //blueRecon.write("C:/Users/ACER/Desktop/blueRecon.png");



    //
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    //


    cout << "\nDeclaring 2nd set of coefficients (6+6+6) ..." << endl;

    // the set of coefficients ai
    float R2_a0 = 0.654108;
    float R2_a1 = -0.0980758;
    float R2_a2 = -0.00294018;
    float R2_a3 = 0.00305956;
    float R2_a4 = 0.44103;
    float R2_a5 = 0.00296609;

    float G2_a0 = 0.652387;
    float G2_a1 = -0.0979924;
    float G2_a2 = -0.00218615;
    float G2_a3 = 0.00133557;
    float G2_a4 = 0.443554;
    float G2_a5 = 0.00285239;

    float B2_a0 = 0.657707;
    float B2_a1 = -0.0998597;
    float B2_a2 = -0.00204011;
    float B2_a3 = -0.00155064;
    float B2_a4 = 0.441307;
    float B2_a5 = 0.00464425;

    cout << "Generating 2nd batch of reconstructed pixels" << endl;

    // testing recon 2
    // 3rd stage of building the recon image, filling in the 2nd batch reconstructed pixels

    for (int i = 1; i <= height2 - 3; i += 2) {                 // i and j is the position of the interested pixel in canva.png
        for (int j = 2; j <= width2 - 4; j += 2) {      // 2 to -4 instead of 1 to -3

            int Ii = ((i * width2) + j) * 3;

            int I0 = (((i - 1) * width2) + j) * 3;
            int I1 = (((i - 1) * width2) + j + 2) * 3;
            int I2 = (((i + 1) * width2) + j) * 3;
            int I3 = (((i + 1) * width2) + j + 2) * 3;
            int I4 = (((i - 1) * width2) + j - 2) * 3;
            int I5 = (((i + 1) * width2) + j - 2) * 3;

            *(Qu2 + Ii) = (*(Qu2 + I0)) * R2_a0 + (*(Qu2 + I1)) * R2_a1 + (*(Qu2 + I2)) * R2_a2 + (*(Qu2 + I3)) * R2_a3 + (*(Qu2 + I4)) * R2_a4 + (*(Qu2 + I5)) * R2_a5;
            *(Qu2 + Ii + 1) = (*(Qu2 + I0 + 1)) * G2_a0 + (*(Qu2 + I1 + 1)) * G2_a1 + (*(Qu2 + I2 + 1)) * G2_a2 + (*(Qu2 + I3 + 1)) * G2_a3 + (*(Qu2 + I4 + 1)) * G2_a4 + (*(Qu2 + I5 + 1)) * G2_a5;
            *(Qu2 + Ii + 2) = (*(Qu2 + I0 + 2)) * B2_a0 + (*(Qu2 + I1 + 2)) * B2_a1 + (*(Qu2 + I2 + 2)) * B2_a2 + (*(Qu2 + I3 + 2)) * B2_a3 + (*(Qu2 + I4 + 2)) * B2_a4 + (*(Qu2 + I5 + 2)) * B2_a5;

            *(Qur + Ii) = (*(Qu2 + I0)) * R2_a0 + (*(Qu2 + I1)) * R2_a1 + (*(Qu2 + I2)) * R2_a2 + (*(Qu2 + I3)) * R2_a3 + (*(Qu2 + I4)) * R2_a4 + (*(Qu2 + I5)) * R2_a5;
            *(Qug + Ii + 1) = (*(Qu2 + I0 + 1)) * G2_a0 + (*(Qu2 + I1 + 1)) * G2_a1 + (*(Qu2 + I2 + 1)) * G2_a2 + (*(Qu2 + I3 + 1)) * G2_a3 + (*(Qu2 + I4 + 1)) * G2_a4 + (*(Qu2 + I5 + 1)) * G2_a5;
            *(Qub + Ii + 2) = (*(Qu2 + I0 + 2)) * B2_a0 + (*(Qu2 + I1 + 2)) * B2_a1 + (*(Qu2 + I2 + 2)) * B2_a2 + (*(Qu2 + I3 + 2)) * B2_a3 + (*(Qu2 + I4 + 2)) * B2_a4 + (*(Qu2 + I5 + 2)) * B2_a5;

        }
        cout << "2nd :  " << i << endl;
    }

    cout << "Done generating the 2nd batch of reconstructed pixels" << endl;

    cout << "Syncing..." << endl;
    canvaPix.sync();
    cout << "Writing image..." << endl;
    canva.write("C:/Users/ACER/Desktop/canva2.png");

    cout << "Syncing RGB..." << endl;
    //redRecPix.sync();
    //greenRecPix.sync();
    //blueRecPix.sync();
    //redRecon.write("C:/Users/ACER/Desktop/redRecon1.png");
    //greenRecon.write("C:/Users/ACER/Desktop/greenRecon1.png");
    //blueRecon.write("C:/Users/ACER/Desktop/blueRecon1.png");

    //
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // so far correct




    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    //



    cout << "\nDeclaring 3rd set of coefficients (6+6+6) ..." << endl;

    // the set of coefficients ai
    float R3_a0 = 0.659152;
    float R3_a1 = -0.00143253;
    float R3_a2 = -0.0990872;
    float R3_a3 = -0.00341;
    float R3_a4 = 0.437125;
    float R3_a5 = 0.00766668;

    float G3_a0 = 0.657795;
    float G3_a1 = -0.0010009;
    float G3_a2 = -0.0986959;
    float G3_a3 = -0.00394662;
    float G3_a4 = 0.438471;
    float G3_a5 = 0.00737687;

    float B3_a0 = 0.66211;
    float B3_a1 = -0.000430273;
    float B3_a2 = -0.0996888;
    float B3_a3 = -0.00529079;
    float B3_a4 = 0.435704;
    float B3_a5 = 0.00758895;



    cout << "Generating 3rd batch of reconstructed pixels" << endl;

    // testing recon 3
    // 4th stage of building the recon image, filling in the 3rd batch reconstructed pixels

    for (int i = 2; i <= height2 - 4; i += 2) {         // 2 to -4 instead of 1 to -3
        for (int j = 1; j <= width2 - 3; j += 2) {          // i and j is the position of the interested pixel in canva.png

            int Ii = ((i * width2) + j) * 3;

            int I0 = ((i * width2) + j - 1) * 3;
            int I1 = ((i * width2) + j + 1) * 3;
            int I2 = (((i + 2) * width2) + j - 1) * 3;
            int I3 = (((i + 2) * width2) + j + 1) * 3;
            int I4 = (((i - 2) * width2) + j - 1) * 3;
            int I5 = (((i - 2) * width2) + j + 1) * 3;

            *(Qu2 + Ii) = (*(Qu2 + I0)) * R3_a0 + (*(Qu2 + I1)) * R3_a1 + (*(Qu2 + I2)) * R3_a2 + (*(Qu2 + I3)) * R3_a3 + (*(Qu2 + I4)) * R3_a4 + (*(Qu2 + I5)) * R3_a5;
            *(Qu2 + Ii + 1) = (*(Qu2 + I0 + 1)) * G3_a0 + (*(Qu2 + I1 + 1)) * G3_a1 + (*(Qu2 + I2 + 1)) * G3_a2 + (*(Qu2 + I3 + 1)) * G3_a3 + (*(Qu2 + I4 + 1)) * G3_a4 + (*(Qu2 + I5 + 1)) * G3_a5;
            *(Qu2 + Ii + 2) = (*(Qu2 + I0 + 2)) * B3_a0 + (*(Qu2 + I1 + 2)) * B3_a1 + (*(Qu2 + I2 + 2)) * B3_a2 + (*(Qu2 + I3 + 2)) * B3_a3 + (*(Qu2 + I4 + 2)) * B3_a4 + (*(Qu2 + I5 + 2)) * B3_a5;

            *(Qur + Ii) = (*(Qu2 + I0)) * R3_a0 + (*(Qu2 + I1)) * R3_a1 + (*(Qu2 + I2)) * R3_a2 + (*(Qu2 + I3)) * R3_a3 + (*(Qu2 + I4)) * R3_a4 + (*(Qu2 + I5)) * R3_a5;
            *(Qug + Ii + 1) = (*(Qu2 + I0 + 1)) * G3_a0 + (*(Qu2 + I1 + 1)) * G3_a1 + (*(Qu2 + I2 + 1)) * G3_a2 + (*(Qu2 + I3 + 1)) * G3_a3 + (*(Qu2 + I4 + 1)) * G3_a4 + (*(Qu2 + I5 + 1)) * G3_a5;
            *(Qub + Ii + 2) = (*(Qu2 + I0 + 2)) * B3_a0 + (*(Qu2 + I1 + 2)) * B3_a1 + (*(Qu2 + I2 + 2)) * B3_a2 + (*(Qu2 + I3 + 2)) * B3_a3 + (*(Qu2 + I4 + 2)) * B3_a4 + (*(Qu2 + I5 + 2)) * B3_a5;
        }
        cout << "3rd :  " << i << endl;
    }

    cout << "Done generating the 3rd batch of reconstructed pixels" << endl;

    cout << "Syncing..." << endl;
    canvaPix.sync();
    cout << "Writing image..." << endl;
    canva.write("C:/Users/ACER/Desktop/canva3.png");

    //cout << "Syncing RGB..." << endl;
    //redRecPix.sync();
    //greenRecPix.sync();
    //blueRecPix.sync();

    //cout << "Writing R..." << endl;
    //redRecon.write("C:/Users/ACER/Desktop/redRecon.png");
    //cout << "Writing G..." << endl;
    //greenRecon.write("C:/Users/ACER/Desktop/greenRecon.png");
    //cout << "Writing B..." << endl;
    //blueRecon.write("C:/Users/ACER/Desktop/blueRecon.png");





    //
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    // 
    //


    return 0;
}



*/










