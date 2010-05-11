#ifndef _ORSA_OSG_FIND_NAMED_NODE_VISITOR_H_
#define _ORSA_OSG_FIND_NAMED_NODE_VISITOR_H_

#include <osg/NodeVisitor>

#include <string>
#include <vector>

namespace orsaOSG {
  
    class FindNamedNodeVisitor : public osg::NodeVisitor {
    
    public:
        FindNamedNodeVisitor(const std::string & name) :
            osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN),
            _name(name) { }
      
    public:
        virtual void apply(osg::Node & node) {
            if (node.getName()==_name) {
                _foundNodes.push_back(&node);
            }
            traverse(node);
        }
    
    public:
        typedef std::vector< osg::ref_ptr<osg::Node> > NodeList;
    
    public:
        std::string _name;
        NodeList _foundNodes;
    };
  
} // namespace orsaOSG

#endif // _ORSA_OSG_FIND_NAMED_NODE_VISITOR_H_
