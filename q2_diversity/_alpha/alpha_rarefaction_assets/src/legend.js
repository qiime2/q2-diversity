import { curData, toggle, toggleAll } from './data';

export default function appendLegendKey(legend, entry, ly, c, color) {
  // line toggle in the legend
  legend.append('rect')
      .attr('id', `idrect${entry}`)
      .attr('class', 'legend')
      .attr('x', 0)
      .attr('y', ly - 2.5)
      .attr('width', 15)
      .attr('height', 5)
      .style('stroke', 'darkGrey')
      .style('fill', () => curData[entry].line)
      .on('click', () => {
        const newColor = curData[entry].line === c ? 'white' : c;
        if (entry === 'Select%20All') {
          toggleAll(false, true, color);
        } else {
          toggle(entry, null, newColor);
          this.style('fill', () => curData[entry].line);
        }
      });
  // dot toggle in the legend
  legend.append('circle')
      .attr('id', `iddot${entry}`)
      .attr('class', 'legend')
      .attr('cx', 30)
      .attr('cy', ly)
      .attr('r', 5)
      .style('stroke', 'darkGrey')
      .style('fill', () => curData[entry].dots)
      .on('click', () => {
        const newColor = curData[entry].dots === c ? 'white' : c;
        if (entry === 'Select%20All') {
          toggleAll(true, false, color);
        } else {
          toggle(entry, newColor, null);
          this.style('fill', () => curData[entry].dots);
        }
      });
  // text for key in the legend
  legend.append('text')
      .attr('class', 'legend')
      .attr('x', 40)
      .attr('y', ly + 5)
      .style('font', '10px sans-serif')
      .text(decodeURI(entry));
}

