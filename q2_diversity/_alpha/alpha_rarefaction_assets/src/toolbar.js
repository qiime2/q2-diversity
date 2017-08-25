/* global d */
import state from './state';

export function addMetricPicker(row, metrics, selectedMetric) {
  row.select('.metricPicker select')
    .on('change', function changeCategory() {
      const newMetric = metrics[this.selectedIndex];
      state.setMetric(newMetric);
    })
    .selectAll('option')
    .data(metrics)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedMetric));
}

export function addCategoryPicker(row, categories, selectedCategory) {
  row.select('.categoryPicker select')
    .on('change', function changeCategory() {
      const newCategory = categories[this.selectedIndex];
      state.setCategory(newCategory);
    })
    .selectAll('option')
    .data(categories)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedCategory));
}
