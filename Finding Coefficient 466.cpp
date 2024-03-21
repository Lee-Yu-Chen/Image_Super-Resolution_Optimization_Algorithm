
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
// 
// 5 Jan 2023
// 
// Agenda : 
// - Explanation of method of implementation of finding coefficients and image reconstruction
// - Demonstration of result, i.e. image reconstruction
// - Review of the quality of reconstructed image and values of coefficients
// 
// Next things to do : 
// - Use a more advance method of implementing the upsampling and finding coeffcient : 
//      - instead of upsampling by 2, umsample by 10
//      - each of the intermediate pixels take references from 9 surrounding pixels from low-res image
// 
// - Start MISR
// - Finding the misalignment between 2 images by fitting a 2-dimensional parabolic curve to each image, 
// - Use what called "Design of Experiments" method to find the coefficients of the parabolic curves
// - then find the distance between the peaks of the two 2-D parabolic curves, and this is where the 2 images correlate the most,
// - then x and y offset can be obtained, then the 2 images can be aligned  
// 
// 
// Notes : 
// - thesis need to include Math
// - Pro forma
// - Image rotation is too complicated for this, maybe suitable for a PhD
// - Correlation should be non-integer
// - It is impossible to get image better than the reference image (the original hi-res image) unless we're using MISR
// - Obtaining significant improvement from the low-res image is already very good 
// - Just keep going and we'll see how far we can go forward, 
// - until one day we'll need to stop and just converge what we've got and write them into the thesis
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



    Image image;
    image.read("C:/Users/ACER/Downloads/Picture/Batu Ferringhi.jpg");
    //image.read("C:/Users/ACER/Desktop/canva3.png");


    size_t w = image.size().width();
    size_t h = image.size().height();
    //image.sample(Geometry(w/2, h/2));
    image.sample(Geometry(w, h));

    image.write("C:/Users/ACER/Desktop/asd.png");




    /*
    // create a 16x16 white image with a red dot
    Image image(Geometry(16, 16), Color("#000000"));         // #7df9ff
    Pixels p(image);
    Quantum* q = p.get(3, 3, 1, 1);

    *q = 65535;
    *(q + 1) = 0;
    *(q + 2) = 0;

    p.sync();
    image.write("C:/Users/ACER/Desktop/sixteen.png");
    *




    w = image.size().width();
    h = image.size().height();


    int s = 2;      // scaling factor

    //Image resized = image;
    //resized.resize(Geometry(width / s, height / s));
    //resized.write("C:/Users/ACER/Desktop/resized.png");

    Image sampled = image;
    sampled.sample(Geometry(w / s, h / s));
    sampled.write("C:/Users/ACER/Desktop/sampled.png");

    //Image scaled = image;
    //scaled.scale(Geometry(width / s, height / s));
    //scaled.write("C:/Users/ACER/Desktop/scaled.png");


    // dimensions of the down-sampled image
    size_t width = sampled.size().width();
    size_t height = sampled.size().height();

    size_t width2 = width * s;
    size_t height2 = height * s;



    // move the existing pixels on the down-sampled image to the reconstruction canva
    Image canva(Geometry(width2, height2), Color("#000000"));     // make an image with the original size
    Pixels imagePix(image);
    Pixels samplePix(sampled);
    Pixels canvaPix(canva);
    Quantum* Qui = imagePix.get(0, 0, w, h);
    Quantum* Qu1 = samplePix.get(0, 0, width, height);
    Quantum* Qu2 = canvaPix.get(0, 0, width2, height2);


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
    canvaPix.sync();
    canva.write("C:/Users/ACER/Desktop/canva.png");





    int* redArr = new int[width * height];
    int* greenArr = new int[width * height];
    int* blueArr = new int[width * height];


    // scale from 65535 (uint16) to 255 (uint8) and store them in dynamic memory array
    // can skip this step if a Matrix of long long int is used
    int I = 0;
    for (int i = 0; i <= height - 1; i++) {
        for (int j = 0; j <= width - 1; j++) {
            *(redArr + I) = (*(Qu1 + I * 3)) * 255 / 65535;
            *(greenArr + I) = (*(Qu1 + I * 3 + 1)) * 255 / 65535;
            *(blueArr + I) = (*(Qu1 + I * 3 + 2)) * 255 / 65535;
            I++;
        }
    }

    // R : 00000045C68FFCE0        00000045C68FFD20        00000045C68FFD60

    // P : 00000045C68FFC60        00000045C68FFC78        00000045C68FFC90



    int* Arr[3] = { redArr, greenArr, blueArr };

    //Eigen::Matrix4f R[3];
    //Eigen::MatrixXf P_R(4, 1), P_G(4, 1), P_B(4, 1);
    //Eigen::MatrixXf P[3] = { P_R, P_G, P_B };

    //Eigen::Matrix<long long int, 4, 4> R[3];
    //Eigen::Matrix<long long int, 4, 1> P_R(4, 1), P_G(4, 1), P_B(4, 1);
    //Eigen::Matrix<long long int, 4, 1> P[3] = { P_R, P_G, P_B };


    Eigen::Matrix<double, 4, 4> R[3];       // R matrices (RGB) for the 1st batch of reconstructed pixels (4x4)
    Eigen::Matrix<double, 6, 6> R2[3];      // R matrices (RGB) for the 2nd batch of reconstructed pixels (6x6)
    Eigen::Matrix<double, 6, 6> R3[3];      // R matrices (RGB) for the 3rd batch of reconstructed pixels (6x6)

    Eigen::Matrix<double, 4, 1> P_R(4, 1), P_G(4, 1), P_B(4, 1);        // P matrices (RGB) for the 1st batch of reconstructed pixels (4x1)
    Eigen::Matrix<double, 6, 1> P_R2(6, 1), P_G2(6, 1), P_B2(6, 1);     // P matrices (RGB) for the 2nd batch of reconstructed pixels (6x1)
    Eigen::Matrix<double, 6, 1> P_R3(6, 1), P_G3(6, 1), P_B3(6, 1);     // P matrices (RGB) for the 3rd batch of reconstructed pixels (6x1)

    Eigen::Matrix<double, 4, 1> P[3] = { P_R, P_G, P_B };       // P matrices array for 3 batches of reconstructed pixels
    Eigen::Matrix<double, 6, 1> P2[3] = { P_R2, P_G2, P_B2 };   // (6x1)
    Eigen::Matrix<double, 6, 1> P3[3] = { P_R3, P_G3, P_B3 };   // (6x1)
    // using float instead of int because only float can perform inverse
    // using double instead of float because memory problem


    cout << "RGBArr address : " << endl;
    cout << redArr << '\t' << greenArr << '\t' << blueArr << endl << endl;

    cout << "Arr[C] address : " << endl;
    cout << Arr[0] << '\t' << Arr[1] << '\t' << Arr[2] << endl << endl;

    cout << "R[C] address : " << endl;
    cout << &(R[0]) << '\t' << &(R[1]) << '\t' << &(R[2]) << endl << endl;

    cout << "P_RGB address : " << endl;
    cout << &(P_R) << '\t' << &(P_G) << '\t' << &(P_B) << endl << endl;

    cout << "P[C] address : " << endl;
    cout << &(P[0]) << '\t' << &(P[1]) << '\t' << &(P[2]) << endl << endl;





    // Reconstruct the 1st set of intermediate pixels, for all 3 channels (RGB)
    // Bottom right pixel
    for (int C = 0; C <= 2; C++) {

        // building the R matrix
        int x0, x1, x2, x3;

        long long int x0x0 = 0, x0x1 = 0, x0x2 = 0, x0x3 = 0;
        long long int x1x0 = 0, x1x1 = 0, x1x2 = 0, x1x3 = 0;
        long long int x2x0 = 0, x2x1 = 0, x2x2 = 0, x2x3 = 0;
        long long int x3x0 = 0, x3x1 = 0, x3x2 = 0, x3x3 = 0;

        for (int i = 0; i <= height - 2; i++) {
            for (int j = 0; j <= width - 2; j++) {              // i and j is the position of the top left pixel in the sampled.png

                x0 = *(Arr[C] + i * width + j);
                //x1 = *(Arr[C] + (i + 1) * width + j);       //
                //x2 = *(Arr[C] + i * width + j + 1);
                x1 = *(Arr[C] + i * width + j + 1);
                x2 = *(Arr[C] + (i + 1) * width + j);       //
                x3 = *(Arr[C] + (i + 1) * width + j + 1);       //

                x0x0 += x0 * x0;
                x0x1 += x0 * x1;
                x0x2 += x0 * x2;
                x0x3 += x0 * x3;

                x1x0 += x1 * x0;
                x1x1 += x1 * x1;
                x1x2 += x1 * x2;
                x1x3 += x1 * x3;

                x2x0 += x2 * x0;
                x2x1 += x2 * x1;
                x2x2 += x2 * x2;
                x2x3 += x2 * x3;

                x3x0 += x3 * x0;
                x3x1 += x3 * x1;
                x3x2 += x3 * x2;
                x3x3 += x3 * x3;

            }
            //cout << "1st C : " << C << "    xi * xi_T     " << i << endl;
        }


        R[C] << x0x0, x0x1, x0x2, x0x3, x1x0, x1x1, x1x2, x1x3, x2x0, x2x1, x2x2, x2x3, x3x0, x3x1, x3x2, x3x3;
        //cout << "1st C : " << C << "   R " << endl;
        //cout << "Address verification : " << &(R[C]) << endl;
        //cout << R[C] << endl << endl;





        // building the P matrix
        int yi;
        long long int x0yi = 0, x1yi = 0, x2yi = 0, x3yi = 0;       // need to change Eigen::MatrixXi to long long int

        for (int i = 0; i <= height - 2; i++) {
            for (int j = 0; j <= width - 2; j++) {          // i and j is the position of the top left pixel in the sampled.png

                x0 = *(Arr[C] + i * width + j);
                //x1 = *(Arr[C] + (i + 1) * width + j);       //
                //x2 = *(Arr[C] + i * width + j + 1);
                x1 = *(Arr[C] + i * width + j + 1);
                x2 = *(Arr[C] + (i + 1) * width + j);
                x3 = *(Arr[C] + (i + 1) * width + j + 1);       //


                int I2 = (((i * 2 + 1) * w) + j * 2 + 1) * 3;               // j*2 no need +1 (?)  // nope, should be fine (j need +1)
                // position of the corresponding pixel in hi-res image with respect to the top left pixel
                // this line is kinda not accurate because w is not equal to width2

                //yi = (*(Qui + I2)) * 255 / 65535;       // + C
                yi = (*(Qui + I2 + C)) * 255 / 65535;

                x0yi += x0 * yi;
                x1yi += x1 * yi;
                x2yi += x2 * yi;
                x3yi += x3 * yi;
            }
            //cout << "1st C : " << C << "    xi * yi     " << i << '\t' << x0yi << '\t' << x1yi << '\t' << x2yi << '\t' << x3yi << endl;
        }


        P[C] << x0yi, x1yi, x2yi, x3yi;
        //cout << "1st C : " << C << "   P " << endl;
        //cout << "Address verification : " << &(P[C]) << endl;
        //cout << P[C] << endl << endl;
    }



    cout << endl << endl << endl;
    for (int i = 0; i <= 2; i++) {
        cout << R[i] << "\n\n" << P[i] << "\n\n\n\n";
    }




    //Eigen::MatrixX4f Red_a(4, 1);
    //Eigen::MatrixX4f Green_a(4,1);
    //Eigen::MatrixX4f Blue_a(4, 1);

    Eigen::Matrix<double, 4, 1> Red_a(4, 1);
    Eigen::Matrix<double, 4, 1> Green_a(4, 1);
    Eigen::Matrix<double, 4, 1> Blue_a(4, 1);


    Red_a = R[0].inverse() * P[0];
    Green_a = R[1].inverse() * P[1];
    Blue_a = R[2].inverse() * P[2];

    cout << "Red_a \n" << Red_a << endl << endl;
    cout << "Green_a \n" << Green_a << endl << endl;
    cout << "Blue_a \n" << Blue_a << endl << endl;




    // the set of coefficients ai
    float R_a0 = Red_a(0, 0);
    float R_a1 = Red_a(1, 0);
    float R_a2 = Red_a(2, 0);
    float R_a3 = Red_a(3, 0);

    float G_a0 = Green_a(0, 0);
    float G_a1 = Green_a(1, 0);
    float G_a2 = Green_a(2, 0);
    float G_a3 = Green_a(3, 0);

    float B_a0 = Blue_a(0, 0);
    float B_a1 = Blue_a(1, 0);
    float B_a2 = Blue_a(2, 0);
    float B_a3 = Blue_a(3, 0);





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
    }

    canvaPix.sync();
    canva.write("C:/Users/ACER/Desktop/canva1.png");

    redRecPix.sync();
    greenRecPix.sync();
    blueRecPix.sync();
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


    // Reconstruct the 2nd set of intermediate pixels, for all 3 channels (RGB)
    // Bottom left pixel
    for (int C = 0; C <= 2; C++) {

        // building the R matrix
        int x0, x1, x2, x3, x4, x5;

        long long int x0x0 = 0, x0x1 = 0, x0x2 = 0, x0x3 = 0, x0x4 = 0, x0x5 = 0;
        long long int x1x0 = 0, x1x1 = 0, x1x2 = 0, x1x3 = 0, x1x4 = 0, x1x5 = 0;
        long long int x2x0 = 0, x2x1 = 0, x2x2 = 0, x2x3 = 0, x2x4 = 0, x2x5 = 0;
        long long int x3x0 = 0, x3x1 = 0, x3x2 = 0, x3x3 = 0, x3x4 = 0, x3x5 = 0;
        long long int x4x0 = 0, x4x1 = 0, x4x2 = 0, x4x3 = 0, x4x4 = 0, x4x5 = 0;
        long long int x5x0 = 0, x5x1 = 0, x5x2 = 0, x5x3 = 0, x5x4 = 0, x5x5 = 0;

        for (int i = 0; i <= height - 2; i++) {                 // i and j is the position of the top left pixel in the sampled.png
            for (int j = 1; j <= width - 2; j++) {          // start from 1 instead of 0

                x0 = *(Arr[C] + i * width + j);
                x1 = *(Arr[C] + i * width + j + 1);
                x2 = *(Arr[C] + (i + 1) * width + j);
                x3 = *(Arr[C] + (i + 1) * width + j + 1);
                x4 = *(Arr[C] + i * width + j - 1);
                x5 = *(Arr[C] + (i + 1) * width + j - 1);

                x0x0 += x0 * x0;
                x0x1 += x0 * x1;
                x0x2 += x0 * x2;
                x0x3 += x0 * x3;
                x0x4 += x0 * x4;
                x0x5 += x0 * x5;

                x1x0 += x1 * x0;
                x1x1 += x1 * x1;
                x1x2 += x1 * x2;
                x1x3 += x1 * x3;
                x1x4 += x1 * x4;
                x1x5 += x1 * x5;

                x2x0 += x2 * x0;
                x2x1 += x2 * x1;
                x2x2 += x2 * x2;
                x2x3 += x2 * x3;
                x2x4 += x2 * x4;
                x2x5 += x2 * x5;

                x3x0 += x3 * x0;
                x3x1 += x3 * x1;
                x3x2 += x3 * x2;
                x3x3 += x3 * x3;
                x3x4 += x3 * x4;
                x3x5 += x3 * x5;

                x4x0 += x4 * x0;
                x4x1 += x4 * x1;
                x4x2 += x4 * x2;
                x4x3 += x4 * x3;
                x4x4 += x4 * x4;
                x4x5 += x4 * x5;

                x5x0 += x5 * x0;
                x5x1 += x5 * x1;
                x5x2 += x5 * x2;
                x5x3 += x5 * x3;
                x5x4 += x5 * x4;
                x5x5 += x5 * x5;
            }
            //cout << "2nd C : " << C << "    xi * xi_T     " << i << endl;
        }


        R2[C] << x0x0, x0x1, x0x2, x0x3, x0x4, x0x5,
            x1x0, x1x1, x1x2, x1x3, x1x4, x1x5,
            x2x0, x2x1, x2x2, x2x3, x2x4, x2x5,
            x3x0, x3x1, x3x2, x3x3, x3x4, x3x5,
            x4x0, x4x1, x4x2, x4x3, x4x4, x4x5,
            x5x0, x5x1, x5x2, x5x3, x5x4, x5x5;

        //cout << "2nd C : " << C << "   R2 " << endl;
        //cout << "Address verification : " << &(R2[C]) << endl;
        //cout << R2[C] << endl << endl;





        // building the P matrix
        int yi;
        long long int x0yi = 0, x1yi = 0, x2yi = 0, x3yi = 0, x4yi = 0, x5yi = 0;       // need to change Eigen::MatrixXi to long long int

        for (int i = 0; i <= height - 2; i++) {             // i and j is the position of the top left pixel in the sampled.png
            for (int j = 1; j <= width - 2; j++) {      // start from 1 instead of 0

                x0 = *(Arr[C] + i * width + j);
                x1 = *(Arr[C] + i * width + j + 1);
                x2 = *(Arr[C] + (i + 1) * width + j);
                x3 = *(Arr[C] + (i + 1) * width + j + 1);
                x4 = *(Arr[C] + i * width + j - 1);
                x5 = *(Arr[C] + (i + 1) * width + j - 1);


                int I2 = (((i * 2 + 1) * w) + j * 2) * 3;   // pixel from the high-res image
                // position of the reconstructed pixel in hi-res image with respect to the top left pixel
                // this line is kinda not accurate because w is not equal to width2

                //yi = (*(Qui + I2)) * 255 / 65535;       // + C
                yi = (*(Qui + I2 + C)) * 255 / 65535;

                x0yi += x0 * yi;
                x1yi += x1 * yi;
                x2yi += x2 * yi;
                x3yi += x3 * yi;
                x4yi += x4 * yi;
                x5yi += x5 * yi;
            }
            //cout << "2nd C : " << C << "    xi * yi     " << i << '\t' << x0yi << '\t' << x1yi << '\t' << x2yi << '\t' << x3yi << '\t' << x4yi
                //<< '\t' << x5yi << endl;
        }


        P2[C] << x0yi, x1yi, x2yi, x3yi, x4yi, x5yi;
        //cout << "2nd C : " << C << "   P2 " << endl;
        //cout << "Address verification : " << &(P2[C]) << endl;
        //cout << P2[C] << endl << endl;
    }


    // print out the R and P matrix
    cout << endl << endl << endl;
    for (int i = 0; i <= 2; i++) {
        cout << R2[i] << "\n\n" << P2[i] << "\n\n\n\n";
    }


    Eigen::Matrix<double, 6, 1> Red2_a(6, 1);
    Eigen::Matrix<double, 6, 1> Green2_a(6, 1);
    Eigen::Matrix<double, 6, 1> Blue2_a(6, 1);


    Red2_a = R2[0].inverse() * P2[0];
    Green2_a = R2[1].inverse() * P2[1];
    Blue2_a = R2[2].inverse() * P2[2];

    cout << "Red2_a \n" << Red2_a << endl << endl;
    cout << "Green2_a \n" << Green2_a << endl << endl;
    cout << "Blue2_a \n" << Blue2_a << endl << endl;




    // the set of coefficients ai
    float R2_a0 = Red2_a(0, 0);
    float R2_a1 = Red2_a(1, 0);
    float R2_a2 = Red2_a(2, 0);
    float R2_a3 = Red2_a(3, 0);
    float R2_a4 = Red2_a(4, 0);
    float R2_a5 = Red2_a(5, 0);

    float G2_a0 = Green2_a(0, 0);
    float G2_a1 = Green2_a(1, 0);
    float G2_a2 = Green2_a(2, 0);
    float G2_a3 = Green2_a(3, 0);
    float G2_a4 = Green2_a(4, 0);
    float G2_a5 = Green2_a(5, 0);

    float B2_a0 = Blue2_a(0, 0);
    float B2_a1 = Blue2_a(1, 0);
    float B2_a2 = Blue2_a(2, 0);
    float B2_a3 = Blue2_a(3, 0);
    float B2_a4 = Blue2_a(4, 0);
    float B2_a5 = Blue2_a(5, 0);


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
    }

    canvaPix.sync();
    canva.write("C:/Users/ACER/Desktop/canva2.png");

    redRecPix.sync();
    greenRecPix.sync();
    blueRecPix.sync();
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


    // Reconstruct the 3rd set of intermediate pixels, for all 3 channels (RGB)
    // Bottom right pixel
    for (int C = 0; C <= 2; C++) {

        // building the R matrix
        int x0, x1, x2, x3, x4, x5;

        long long int x0x0 = 0, x0x1 = 0, x0x2 = 0, x0x3 = 0, x0x4 = 0, x0x5 = 0;
        long long int x1x0 = 0, x1x1 = 0, x1x2 = 0, x1x3 = 0, x1x4 = 0, x1x5 = 0;
        long long int x2x0 = 0, x2x1 = 0, x2x2 = 0, x2x3 = 0, x2x4 = 0, x2x5 = 0;
        long long int x3x0 = 0, x3x1 = 0, x3x2 = 0, x3x3 = 0, x3x4 = 0, x3x5 = 0;
        long long int x4x0 = 0, x4x1 = 0, x4x2 = 0, x4x3 = 0, x4x4 = 0, x4x5 = 0;
        long long int x5x0 = 0, x5x1 = 0, x5x2 = 0, x5x3 = 0, x5x4 = 0, x5x5 = 0;

        for (int i = 1; i <= height - 2; i++) {               // start from 1 instead of 0
            for (int j = 0; j <= width - 2; j++) {          // i and j is the position of the top left pixel in the sampled.png


                x0 = *(Arr[C] + i * width + j);
                x1 = *(Arr[C] + i * width + j + 1);
                x2 = *(Arr[C] + (i + 1) * width + j);
                x3 = *(Arr[C] + (i + 1) * width + j + 1);
                x4 = *(Arr[C] + (i - 1) * width + j);
                x5 = *(Arr[C] + (i - 1) * width + j + 1);

                x0x0 += x0 * x0;
                x0x1 += x0 * x1;
                x0x2 += x0 * x2;
                x0x3 += x0 * x3;
                x0x4 += x0 * x4;
                x0x5 += x0 * x5;

                x1x0 += x1 * x0;
                x1x1 += x1 * x1;
                x1x2 += x1 * x2;
                x1x3 += x1 * x3;
                x1x4 += x1 * x4;
                x1x5 += x1 * x5;

                x2x0 += x2 * x0;
                x2x1 += x2 * x1;
                x2x2 += x2 * x2;
                x2x3 += x2 * x3;
                x2x4 += x2 * x4;
                x2x5 += x2 * x5;

                x3x0 += x3 * x0;
                x3x1 += x3 * x1;
                x3x2 += x3 * x2;
                x3x3 += x3 * x3;
                x3x4 += x3 * x4;
                x3x5 += x3 * x5;

                x4x0 += x4 * x0;
                x4x1 += x4 * x1;
                x4x2 += x4 * x2;
                x4x3 += x4 * x3;
                x4x4 += x4 * x4;
                x4x5 += x4 * x5;

                x5x0 += x5 * x0;
                x5x1 += x5 * x1;
                x5x2 += x5 * x2;
                x5x3 += x5 * x3;
                x5x4 += x5 * x4;
                x5x5 += x5 * x5;
            }
            //cout << "3rd C : " << C << "    xi * xi_T     " << i << endl;
        }


        R3[C] << x0x0, x0x1, x0x2, x0x3, x0x4, x0x5,
            x1x0, x1x1, x1x2, x1x3, x1x4, x1x5,
            x2x0, x2x1, x2x2, x2x3, x2x4, x2x5,
            x3x0, x3x1, x3x2, x3x3, x3x4, x3x5,
            x4x0, x4x1, x4x2, x4x3, x4x4, x4x5,
            x5x0, x5x1, x5x2, x5x3, x5x4, x5x5;

        //cout << "3rd C : " << C << "   R3 " << endl;
        //cout << "Address verification : " << &(R3[C]) << endl;
        //cout << R3[C] << endl << endl;





        // building the P matrix
        int yi;
        long long int x0yi = 0, x1yi = 0, x2yi = 0, x3yi = 0, x4yi = 0, x5yi = 0;       // need to change Eigen::MatrixXi to long long int

        for (int i = 1; i <= height - 2; i++) {        // start from 1 instead of 0
            for (int j = 0; j <= width - 2; j++) {      // i and j is the position of the top left pixel in the sampled.png

                x0 = *(Arr[C] + i * width + j);
                x1 = *(Arr[C] + i * width + j + 1);
                x2 = *(Arr[C] + (i + 1) * width + j);
                x3 = *(Arr[C] + (i + 1) * width + j + 1);
                x4 = *(Arr[C] + (i - 1) * width + j);
                x5 = *(Arr[C] + (i - 1) * width + j + 1);


                int I2 = (((i * 2) * w) + j * 2 + 1) * 3;   // pixel from the high-res image
                // position of the reconstructed pixel in hi-res image with respect to the top left pixel
                // this line is kinda not accurate because w is not equal to width2

                //yi = (*(Qui + I2)) * 255 / 65535;       // + C
                yi = (*(Qui + I2 + C)) * 255 / 65535;

                x0yi += x0 * yi;
                x1yi += x1 * yi;
                x2yi += x2 * yi;
                x3yi += x3 * yi;
                x4yi += x4 * yi;
                x5yi += x5 * yi;
            }
            //cout << "3rd C : " << C << "    xi * yi     " << i << '\t' << x0yi << '\t' << x1yi << '\t' << x2yi << '\t' << x3yi << '\t' << x4yi
                //<< '\t' << x5yi << endl;
        }


        P3[C] << x0yi, x1yi, x2yi, x3yi, x4yi, x5yi;
        //cout << "3rd C : " << C << "   P3 " << endl;
        //cout << "Address verification : " << &(P3[C]) << endl;
        //cout << P3[C] << endl << endl;
    }


    // print out the R and P matrix
    cout << endl << endl << endl;
    for (int i = 0; i <= 2; i++) {
        cout << R3[i] << "\n\n" << P3[i] << "\n\n\n\n";
    }


    Eigen::Matrix<double, 6, 1> Red3_a(6, 1);
    Eigen::Matrix<double, 6, 1> Green3_a(6, 1);
    Eigen::Matrix<double, 6, 1> Blue3_a(6, 1);


    Red3_a = R3[0].inverse() * P3[0];
    Green3_a = R3[1].inverse() * P3[1];
    Blue3_a = R3[2].inverse() * P3[2];

    cout << "Red3_a \n" << Red3_a << endl << endl;
    cout << "Green3_a \n" << Green3_a << endl << endl;
    cout << "Blue3_a \n" << Blue3_a << endl << endl;



    // the set of coefficients ai
    float R3_a0 = Red3_a(0, 0);
    float R3_a1 = Red3_a(1, 0);
    float R3_a2 = Red3_a(2, 0);
    float R3_a3 = Red3_a(3, 0);
    float R3_a4 = Red3_a(4, 0);
    float R3_a5 = Red3_a(5, 0);

    float G3_a0 = Green3_a(0, 0);
    float G3_a1 = Green3_a(1, 0);
    float G3_a2 = Green3_a(2, 0);
    float G3_a3 = Green3_a(3, 0);
    float G3_a4 = Green3_a(4, 0);
    float G3_a5 = Green3_a(5, 0);

    float B3_a0 = Blue3_a(0, 0);
    float B3_a1 = Blue3_a(1, 0);
    float B3_a2 = Blue3_a(2, 0);
    float B3_a3 = Blue3_a(3, 0);
    float B3_a4 = Blue3_a(4, 0);
    float B3_a5 = Blue3_a(5, 0);


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
    }

    canvaPix.sync();
    canva.write("C:/Users/ACER/Desktop/canva3.png");

    redRecPix.sync();
    greenRecPix.sync();
    blueRecPix.sync();
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



    // compare and compute cost function




    delete redArr;
    delete greenArr;
    delete blueArr;



    return 0;
}

*/














