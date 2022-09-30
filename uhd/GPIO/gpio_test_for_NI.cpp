#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <iostream>
#include <thread>

int UHD_SAFE_MAIN(int argc, char* argv[])
{

    uhd::device_addr_t my_usrp_ip;
    my_usrp_ip["addr0"] = "192.168.30.2";
    uhd::usrp::multi_usrp::sptr  my_usrp = uhd::usrp::multi_usrp::make(my_usrp_ip);
    std::cout << boost::format("Using Device: %s") % my_usrp->get_pp_string() << std::endl;

    std::cout << "Available GPIO banks: " << std::endl;
    auto banks = my_usrp->get_gpio_banks(0);
    for (auto& bank : banks) {
        std::cout << "* " << bank << std::endl;
    }
    
    // set up some catch-all masks
    uint32_t gpio_line = 0xF; // only the bottom 4 lines: 0xF = 00001111 = Pin 0, 1, 2, 3
    uint32_t all_one = 0xFF;
    uint32_t all_zero = 0x00;

    // reset usrp time to 0.00
    my_usrp->set_time_source("internal");
    my_usrp->set_time_next_pps(uhd::time_spec_t(0.0));
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    uhd::time_spec_t now_time = my_usrp->get_time_last_pps();        // define t=0
    // set gpio pins up for output
    my_usrp->set_gpio_attr("FP0", "DDR", all_one, gpio_line, 0);
    my_usrp->set_gpio_attr("FP0", "CTRL", all_zero, gpio_line, 0);
    double Frequency = 42000; //  khz
    double time_step = 1 / Frequency;
    int microseconds_constant = 1000000;
    int sleep_contant_microsec = time_step*microseconds_constant;
    // try{}
    while(true){
        // uhd::time_spec_t now_time = my_usrp->get_time_last_pps();        // define t=0
        for(int i=0;i<Frequency*Frequency;i+=2){ // "Frequency" secs....
            // set all gpio lines to output 1
            my_usrp->clear_command_time();
            my_usrp->set_command_time(now_time + uhd::time_spec_t(time_step*i));
            my_usrp->set_gpio_attr("FP0", "OUT", all_one, gpio_line, 0);     // reset HIGH (async)
            // set all gpio lines to output 0
            my_usrp->clear_command_time();
            my_usrp->set_command_time(now_time + uhd::time_spec_t(time_step*i+time_step));
            my_usrp->set_gpio_attr("FP0", "OUT", all_zero, gpio_line, 0);    // set LOW @ t=2
            // passive waiting
            // std::this_thread::sleep_for(std::chrono::microseconds(sleep_contant_microsec));
            // active waiting
            auto end = std::chrono::steady_clock::now() + std::chrono::microseconds(sleep_contant_microsec);
            while(std::chrono::steady_clock::now() < end);
        }   
        std::cout<<"reach end..."<<std::endl;
        
    }

    return EXIT_SUCCESS;
}

