//
// Created by zliu58 on 6/27/17.
//

#include<fstream>
#include <algorithm>
#include <string>
#include <set>
#include "RouteAnalysis.h"

#include <tulip/GlScene.h>
#include <tulip/BooleanProperty.h>
#include <tulip/StringProperty.h>
#include <tulip/IntegerProperty.h>
#include <tulip/ColorProperty.h>
#include "fabric.h"
#include "ibautils/ib_fabric.h"
#include "ibautils/ib_parser.h"
#include "ibautils/ib_port.h"
#include "ibautils/regex.h"

using namespace std;
using namespace tlp;

PLUGIN(RouteAnalysis)

static const char * paramHelp[] = {
        // File to Open
        HTML_HELP_OPEN() \
  HTML_HELP_DEF( "type", "pathname" ) \
  HTML_HELP_BODY() \
  "Path to ibdiagnet2.fdbs file to import" \
  HTML_HELP_CLOSE()
};

RouteAnalysis::RouteAnalysis(tlp::PluginContext* context)
        : tlp::Algorithm(context)
{
    addInParameter<std::string>("file::filename", paramHelp[0],"");
}


namespace ib = infiniband;
namespace ibp = infiniband::parser;

bool RouteAnalysis::run(){
    assert(graph);

    static const size_t STEPS = 5;
    if(pluginProgress)
    {
        pluginProgress->showPreview(false);
        pluginProgress->setComment("Starting to Import Routes");
        pluginProgress->progress(0, STEPS);
    }

    /**
     * while this does not import
     * nodes/edges, it imports properties
     * for an existing fabric
     */

    ib::tulip_fabric_t * const fabric = ib::tulip_fabric_t::find_fabric(graph, false);
    if(!fabric)
    {
        if(pluginProgress)
            pluginProgress->setError("Unable find fabric. Make sure to preserve data when importing data.");

        return false;
    }

    if(pluginProgress)
    {
        pluginProgress->setComment("Found Fabric");
        pluginProgress->progress(1, STEPS);
    }

    /**
     * Open file to read and import per type
     */
    std::string filename;

    dataSet->get("file::filename", filename);
    std::ifstream ifs(filename.c_str());
    if(!ifs)
    {
        if(pluginProgress)
            pluginProgress->setError("Unable open source file.");

        return false;
    }

    if(pluginProgress)
    {
        pluginProgress->progress(2, STEPS);
        pluginProgress->setComment("Parsing Routes.");
    }

    ibp::ibdiagnet_fwd_db parser;
    if(!parser.parse(*fabric, ifs))
    {
        if(pluginProgress)
            pluginProgress->setError("Unable parse routes file.");

        return false;
    }

    if(pluginProgress)
    {
        pluginProgress->setComment("Parsing Routes complete.");
        pluginProgress->progress(3, STEPS);
    }

    ifs.close();


    //find the source and target nodes of the path
    BooleanProperty *selectBool = graph->getLocalProperty<BooleanProperty>("viewSelection");
    StringProperty *getGuid = graph->getLocalProperty<StringProperty>("ibGuid");
    IntegerProperty *getPortNum = graph->getLocalProperty<IntegerProperty>("ibPortNum");
    ColorProperty *setColor = graph->getLocalProperty<ColorProperty>("viewColor");
    vector<tlp::node> nodes_guid;

    tlp::Iterator<tlp::node> *selections = selectBool->getNodesEqualTo(true,NULL);
    const ib::fabric_t::entities_t &entities_map = fabric->get_entities();
     
        cout<<"Test1"<<endl;
    while(selections->hasNext()){
            cout<<"Test0"<<endl;
        const tlp::node &mynode = selections->next();
        nodes_guid.push_back(mynode);
    }

        cout<<"Test Vector Size: "<<nodes_guid.size()<<endl;
    //const unsigned long long int key1 = std::stol((getGuid->getNodeStringValue(nodes_guid[0])).c_str(),NULL,0);
    //const unsigned long long int key2 = std::stol((getGuid->getNodeStringValue(nodes_guid[1])).c_str(),NULL,0);

    const ib::entity_t & source_node = entities_map.find(std::stol((getGuid->getNodeStringValue(nodes_guid[0])).c_str(),NULL,0))->second;
    const ib::entity_t & target_node = entities_map.find(std::stol((getGuid->getNodeStringValue(nodes_guid[1])).c_str(),NULL,0))->second;

    cout<<"source_guid: "<<source_node.guid<<endl;
    cout<<"target_guid: "<<target_node.guid<<endl;
    cout<<" -------------------------------- "<<endl;
    /*ib::lid_t target_lid = target_node.lid();
    cout<<"target_lid: "<<target_lid<<endl;
    cout<<"target_guid: "<<target_node.guid<<endl;

    std::vector<ib::entity_t *> tmp;
    tmp.push_back(const_cast<ib::entity_t *> (&source_node));

    const ib::entity_t & temp = *tmp.back();
    cout<<"source_guid: "<<source_node.guid<<endl;
    cout<< "Get into the loop"<<endl;*/
        
    ib::lid_t target_lid = target_node.lid();
    std::vector<ib::entity_t *> tmp;

    // Use to count the number of hops;
    unsigned int count_hops = 0;
    
    // Whether the souce node is HCA or not
    if(getPortNum->getNodeValue(nodes_guid[0])== 1){
        //find the only port in HCA
        const ib::entity_t::portmap_t::const_iterator Myport = source_node.ports.begin();

        //use the typedef std::map<port_t*, tlp::edge> port_edges_t to find the edge
        ib::tulip_fabric_t::port_edges_t::iterator Myedge = fabric->port_edges.find(Myport->second);
        setColor->setEdgeValue(Myedge->second, tlp::Color::SpringGreen);
        selectBool->setEdgeValue(Myedge->second, true);
        if(graph->source(Myedge->second).id == nodes_guid[0].id){
            const tlp::edge &e = Myedge->second;
            tmp.push_back(const_cast<ib::entity_t *> ( & entities_map.find(std::stol(getGuid->getNodeStringValue(graph->target(e)).c_str(),NULL,0))->second));
        }else{
            tmp.push_back(const_cast<ib::entity_t *> (&source_node));
        }
    }else{
   
            tmp.push_back(const_cast<ib::entity_t *> (&source_node));
    }
    
    //Whether the target is HCE or not 
    if(getPortNum->getNodeValue(nodes_guid[1])== 1){
        //find the only port in HCA
        const ib::entity_t::portmap_t::const_iterator Myport = target_node.ports.begin();

        //use the typedef std::map<port_t*, tlp::edge> port_edges_t to find the edge
        ib::tulip_fabric_t::port_edges_t::iterator Myedge = fabric->port_edges.find(Myport->second);
        setColor->setEdgeValue(Myedge->second, tlp::Color::SpringGreen);
        selectBool->setEdgeValue(Myedge->second, true);
        if(graph->source(Myedge->second).id == nodes_guid[0].id){
            const tlp::edge &e = Myedge->second;
            const ib::entity_t & real_target = entities_map.find(std::stol((getGuid->getNodeStringValue(graph->target(e))).c_str(),NULL,0))->second;
            while(tmp.back()->guid!= real_target.guid) {
                const ib::entity_t & temp = *tmp.back();
                cout<<"The "<<count_hops<<" step: "<<temp.guid<<endl;
                for (
                        ib::entity_t::routes_t::const_iterator
                                ritr = temp.get_routes().begin(),
                                reitr = temp.get_routes().end();
                        ritr != reitr;
                        ++ritr
                        ) {
                    std::set<ib::lid_t>::const_iterator itr = ritr->second.find(target_lid);
                    if (itr != ritr->second.end()) {
                        const ib::entity_t::portmap_t::const_iterator port_itr = temp.ports.find(ritr->first);
                        if (port_itr != temp.ports.end()) {
                            const ib::port_t *const port = port_itr->second;
                            const ib::tulip_fabric_t::port_edges_t::const_iterator edge_itr = fabric->port_edges.find(
                                    const_cast<ib::port_t *>(port));
                            if (edge_itr != fabric->port_edges.end()) {
                                const tlp::edge &edge = edge_itr->second;
                                setColor->setEdgeValue(edge, tlp::Color::SpringGreen);
                                selectBool->setEdgeValue(edge, true);
                                const ib::entity_t &node = entities_map.find(std::stol((getGuid->getNodeStringValue(graph->target(edge))).c_str(), NULL, 0))->second;
                                tmp.push_back(const_cast<ib::entity_t *> (&node));
                                    count_hops++;
                                    }
                                }
                            }
                        }
              }
              cout<<"The "<<count_hops<<" step: "<<real_target.guid<<endl;
                
        }else{
            //find the only port in HCA
            const ib::entity_t::portmap_t::const_iterator Myport = target_node.ports.begin();

            //use the typedef std::map<port_t*, tlp::edge> port_edges_t to find the edge
            ib::tulip_fabric_t::port_edges_t::iterator Myedge = fabric->port_edges.find(Myport->second);
            const tlp::edge &e = Myedge->second;
            setColor->setEdgeValue(e, tlp::Color::SpringGreen);
            selectBool->setEdgeValue(e, true);
            const ib::entity_t & real_target = entities_map.find(std::stol((getGuid->getNodeStringValue(graph->source(e))).c_str(),NULL,0))->second;
            while(tmp.back()->guid!= real_target.guid) {
                const ib::entity_t & temp = *tmp.back();
                cout<<"The "<<count_hops<<" step: "<<temp.guid<<endl;
                for (
                        ib::entity_t::routes_t::const_iterator
                                ritr = temp.get_routes().begin(),
                                reitr = temp.get_routes().end();
                        ritr != reitr;
                        ++ritr
                        ) {
                    
                    std::set<ib::lid_t>::const_iterator itr = ritr->second.find(target_lid);
                    if (itr != ritr->second.end()) {
                        const ib::entity_t::portmap_t::const_iterator port_itr = temp.ports.find(ritr->first);
                        if (port_itr != temp.ports.end()) {
                            const ib::port_t *const port = port_itr->second;
                            const ib::tulip_fabric_t::port_edges_t::const_iterator edge_itr = fabric->port_edges.find(
                                    const_cast<ib::port_t *>(port));
                            if (edge_itr != fabric->port_edges.end()) {
                                const tlp::edge &edge = edge_itr->second;
                                setColor->setEdgeValue(edge, tlp::Color::SpringGreen);
                                selectBool->setEdgeValue(edge, true);
                                const ib::entity_t &node = entities_map.find(std::stol((getGuid->getNodeStringValue(graph->target(edge))).c_str(), NULL, 0))->second;
                                tmp.push_back(const_cast<ib::entity_t *> (&node));
                                count_hops++;
                            }
                         }
                     }else{
                    
                            cout<<"Do not find routesmap"<<endl;
                    }
                 }
              } 
              cout<<"The "<<count_hops<<" step: "<<real_target.guid<<endl;
                
        }
            
            
        cout<<"The "<<++count_hops<<" step: "<<target_node.guid<<endl; 
        cout<<" ------------------------" <<endl;
        cout<<"The total hops: "<<count_hops<<endl;
    }
    else{
            
        while(tmp.back()->guid!= target_node.guid) {
        const ib::entity_t *temp = tmp.back();
        cout<<"The "<<count_hops<<" step: "<<temp->guid<<endl;
        for (
                ib::entity_t::routes_t::const_iterator
                        ritr = temp->get_routes().begin(),
                        reitr = temp->get_routes().end();
                ritr != reitr;
                ++ritr
                ) {
            std::set<ib::lid_t>::const_iterator itr = ritr->second.find(target_lid);
            if (itr != ritr->second.end()) {
                const ib::entity_t::portmap_t::const_iterator port_itr = temp->ports.find(ritr->first);
                if (port_itr != temp->ports.end()) {
                    const ib::port_t *const port = port_itr->second;
                    const ib::tulip_fabric_t::port_edges_t::const_iterator edge_itr = fabric->port_edges.find(
                            const_cast<ib::port_t *>(port));
                    if (edge_itr != fabric->port_edges.end()) {
                        const tlp::edge &edge = edge_itr->second;
                        //setColor->setEdgeValue(edge, tlp::Color::SpringGreen);
                        selectBool->setEdgeValue(edge, true);
                        const ib::entity_t &node = entities_map.find(std::stol((getGuid->getNodeStringValue(graph->target(edge))).c_str(), NULL, 0))->second;
                        tmp.push_back(const_cast<ib::entity_t *> (&node));
                            count_hops++;
                            }
                        }
                    }
                }
            }


    cout<<"The "<<count_hops<<" step: "<<target_node.guid<<endl;
    cout<<" ------------------------" <<endl;
    cout<<"The total hops: "<<count_hops<<endl;
           
    }
        
       

    if(pluginProgress)
    {
        pluginProgress->setComment("Print the Real hops");
        pluginProgress->progress(STEPS, STEPS);
    }

    return true;
}

